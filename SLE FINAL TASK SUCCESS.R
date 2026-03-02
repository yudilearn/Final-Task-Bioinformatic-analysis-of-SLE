# ==============================================================================
# MASTER SCRIPT: Analisis Ekspresi Gen SLE (Systemic Lupus Erythematosus)
# Dataset: GSE20864 | Platform: Agilent GPL1291
# Tujuan: Identifikasi DEG & Interferon Signature
# ==============================================================================

# --- PART 1: LOAD LIBRARIES ---
# Memastikan semua library yang dibutuhkan tersedia
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(AnnotationDbi)
library(umap)

# --- PART 2: PENGAMBILAN & PRE-PROCESSING DATA ---
message("Sedang mengunduh data dari GEO...")
gset <- getGEO("GSE20864", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

# Matriks ekspresi
ex <- exprs(gset)

# Log2 Transformasi (Hanya jika data belum di-log)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
  message("Log2 transformasi berhasil dilakukan.")
}

# --- PART 3: DEFINISI KELOMPOK SAMPEL ---
metadata <- pData(gset)
healthy_ids <- c("GSM443427", "GSM443428", "GSM443429", "GSM443430", "GSM443431", 
                 "GSM443432", "GSM443433", "GSM443434", "GSM443435", "GSM443436",  
                 "GSM443437", "GSM443438", "GSM443439", "GSM443440", "GSM443441", 
                 "GSM443442", "GSM443443", "GSM443444", "GSM443445", "GSM443446",  
                 "GSM443447", "GSM443448", "GSM443449", "GSM443450", "GSM443451", 
                 "GSM443452", "GSM443453", "GSM443454", "GSM443455", "GSM443456",  
                 "GSM443457", "GSM443458", "GSM443459", "GSM443460", "GSM443461", 
                 "GSM443462", "GSM443463", "GSM443464", "GSM443465", "GSM443466",  
                 "GSM443467", "GSM443468", "GSM443469", "GSM443470", "GSM443471")

# Klasifikasi Healthy vs SLE
metadata$group <- ifelse(metadata$geo_accession %in% healthy_ids, "Healthy", "SLE")
gset$group <- factor(metadata$group, levels = c("Healthy", "SLE"))

#Menjalankan PCA
message("Menjalankan analisis PCA...")

# 1. Transpose matriks (R butuh sampel di baris, gen di kolom)
# Kita gunakan top 500 gen dengan varians tertinggi agar plot lebih bersih
v <- apply(ex, 1, var)
ex_top <- ex[order(v, decreasing = TRUE)[1:500], ]
pca_data <- t(ex_top)

# 2. Hitung PCA
pca_res <- prcomp(pca_data, scale. = TRUE)

# 3. Siapkan data untuk plotting
pca_df <- as.data.frame(pca_res$x)
pca_df$Group <- metadata$group

# 4. Visualisasi dengan ggplot2
library(ggplot2)
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  scale_color_manual(values = c("Healthy" = "blue", "SLE" = "red")) +
  labs(title = "PCA Plot: Healthy vs SLE (GSE20864)",
       subtitle = "Menggunakan Top 500 Gen dengan Varians Tertinggi",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)"))

# Tampilkan dan Simpan
print(p_pca)
ggsave("PCA_GSE20864.png", p_pca, width = 8, height = 6, dpi = 300)

# --- PART 4: ANALISIS STATISTIK (LIMMA) ---
# 1. Design Matrix & Pembersihan Data
design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)

# Imputasi NA & Hapus varians nol
if(any(is.na(ex))) ex <- apply(ex, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))
ex <- ex[apply(ex, 1, var) > 0, ]

# 2. Linear Modeling
fit <- lmFit(ex, design)
contrast_matrix <- makeContrasts(contrasts = "SLE-Healthy", levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 3. Hasil DEG Awal
topTableResults <- topTable(fit2, adjust = "fdr", number = Inf)
topTableResults$PROBEID <- rownames(topTableResults)

# --- PART 5: ANOTASI GEN (AGILENT FIX) ---
feature_info <- fData(gset)
topTableResults$SYMBOL <- feature_info$`Gene symbol`[match(topTableResults$PROBEID, feature_info$ID)]
topTableResults$GENENAME <- feature_info$`Gene title`[match(topTableResults$PROBEID, feature_info$ID)]

# Simpan versi bersih tanpa NA untuk visualisasi
topTable_clean <- topTableResults %>% filter(!is.na(SYMBOL) & SYMBOL != "")

# --- PART 6: VISUALISASI ---

# A. UMAP (Melihat pengelompokan sampel)
umap_out <- umap(t(ex))
df_umap <- data.frame(UMAP1 = umap_out$layout[,1], UMAP2 = umap_out$layout[,2], Group = gset$group)
ggplot(df_umap, aes(UMAP1, UMAP2, color = Group)) + geom_point(size=3) + 
  theme_minimal() + labs(title="Pemisahan Sampel SLE vs Healthy (UMAP)")

# B. Volcano Plot
volcano_df <- topTableResults %>%
  mutate(status = case_when(logFC > 1 & adj.P.Val < 0.05 ~ "UP",
                            logFC < -1 & adj.P.Val < 0.05 ~ "DOWN",
                            TRUE ~ "NO"))

ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.4) + scale_color_manual(values = c("blue", "grey", "red")) +
  theme_minimal() + labs(title="Gen Naik/Turun pada SLE")

# C. Heatmap Top 50 Gen
# Sinkronisasi Matriks agar tidak error "subscript out of bounds"
top50_list <- topTable_clean %>% arrange(adj.P.Val) %>% head(50)
probes_plot <- intersect(top50_list$PROBEID, rownames(ex))
mat_heatmap <- ex[probes_plot, ]

# Penamaan baris yang unik (Gen Symbol)
rownames(mat_heatmap) <- make.unique(as.character(top50_list$SYMBOL[match(probes_plot, top50_list$PROBEID)]))

pheatmap(mat_heatmap, scale = "row", 
         annotation_col = data.frame(Group = gset$group, row.names = colnames(ex)),
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top 50 Gen: Interferon Signature Detected", show_colnames = FALSE)

# --- PART 7: PENYIMPANAN ---
write.csv(topTableResults, "Hasil_Analisis_SLE_GSE20864_Lengkap.csv", row.names = FALSE)
message("Proses Selesai. Hasil telah disimpan dalam CSV.")

# --- PART 8: EKSTRAKSI GEN SIGNIFIKAN (UP & DOWN) ---

# 1. Ekstrak daftar gen yang UP-regulated (Naik pada pasien SLE)
# Menggunakan filter logFC > 1 (kenaikan 2 kali lipat) dan adj.P.Val < 0.01
genes_up <- topTableResults %>%
  filter(logFC > 1 & adj.P.Val < 0.01 & !is.na(SYMBOL)) %>%
  arrange(desc(logFC)) 

# 2. Ekstrak daftar gen yang DOWN-regulated (Turun pada pasien SLE)
# Menggunakan filter logFC < -1 dan adj.P.Val < 0.01
genes_down <- topTableResults %>%
  filter(logFC < -1 & adj.P.Val < 0.01 & !is.na(SYMBOL)) %>%
  arrange(logFC) 

# 3. Menampilkan ringkasan statistik
cat("\n--- RINGKASAN HASIL ANALISIS ---\n")
message("Jumlah gen yang Up-regulated  : ", nrow(genes_up))
message("Jumlah gen yang Down-regulated: ", nrow(genes_down))

# 4. Menampilkan 10 gen teratas untuk masing-masing kategori
print("Top 10 Gen Up-regulated (Potensial Biomarker SLE):")
print(head(genes_up[, c("SYMBOL", "logFC", "adj.P.Val")], 10))

print("Top 10 Gen Down-regulated:")
print(head(genes_down[, c("SYMBOL", "logFC", "adj.P.Val")], 10))

# 5. Simpan ke file CSV terpisah
write.csv(genes_up, "Daftar_Gen_UP_SLE.csv", row.names = FALSE)
write.csv(genes_down, "Daftar_Gen_DOWN_SLE.csv", row.names = FALSE)


# --- PART 9: ANALISIS FUNGSIONAL (GO & KEGG) ---
library(clusterProfiler)
library(org.Hs.eg.db) # Database anotasi manusia
library(enrichplot)


# 1. Konversi Symbol ke Entrez ID dengan Pembersihan Spasi
genes_up_clean <- genes_up %>%
  mutate(SYMBOL = trimws(SYMBOL)) %>% # Menghapus spasi di awal/akhir
  filter(SYMBOL != "")

gen_up_entrez <- bitr(genes_up_clean$SYMBOL, 
                      fromType = "SYMBOL", 
                      toType   = "ENTREZID", 
                      OrgDb    = org.Hs.eg.db)

# Gabungkan logFC ke data Entrez untuk visualisasi cnetplot nanti
gen_up_entrez <- gen_up_entrez %>%
  left_join(genes_up_clean[, c("SYMBOL", "logFC")], by = "SYMBOL")

# 2. Gene Ontology (GO) - Biological Process
# Kita gunakan gen yang berhasil di-map (90% sisanya sudah sangat cukup)
ego_up <- enrichGO(gene          = gen_up_entrez$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP", 
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)

# 3. KEGG Pathway Enrichment
ekegg_up <- enrichKEGG(gene         = gen_up_entrez$ENTREZID,
                       organism     = 'hsa', 
                       pvalueCutoff = 0.05)

# --- PART 10: VISUALISASI HASIL ---

# A. Dotplot GO
# Menampilkan jalur sistem imun yang dominan pada SLE
dotplot(ego_up, showCategory = 15) + 
  labs(title = "GO Enrichment: Jalur Biologis Teraktif pada SLE")
plot_go <- dotplot(ego_up, showCategory = 15) + 
  labs(title = "GO Enrichment: Jalur Biologis Teraktif pada SLE")
#Penyimpanan file
ggsave("7_GO_Enrichment_SLE_Up.png", 
       plot = plot_go, 
       width = 10, 
       height = 8, 
       dpi = 300)

cat("File '7_GO_Enrichment_SLE_Up.png' telah berhasil disimpan.\n")

# B. Dotplot KEGG
# Jika ekegg_up tidak kosong, visualisasikan
if(!is.null(ekegg_up) && nrow(ekegg_up) > 0) {
  dotplot(ekegg_up, showCategory = 10) + labs(title = "KEGG Pathway: Mekanisme Penyakit")
}
plot_kegg <- dotplot(ekegg_up, showCategory = 10) + 
  labs(title = "KEGG Pathway: Mekanisme Penyakit pada SLE")

# Tampilkan di RStudio
print(plot_kegg)

# 2. Simpan ke file gambar
ggsave("8_KEGG_Pathways_SLE_Up.png", 
       plot = plot_kegg, 
       width = 10, 
       height = 8, 
       dpi = 300)
cat("Berhasil! File '8_KEGG_Pathways_SLE_Up.png' telah disimpan.\n")

# C. Cnetplot (Interaksi Gen-Jalur)
# Menghubungkan gen UP-regulated dengan jalur "Response to Type I Interferon"
# Kita ambil logFC agar warna titik gen mencerminkan tingkat kenaikannya
# --- PERBAIKAN PART 10.C: CNETPLOT ---


# 1. Pastikan vektor logFC bersih
gene_list_plot <- gen_up_entrez$logFC
names(gene_list_plot) <- gen_up_entrez$SYMBOL
gene_list_plot <- gene_list_plot[!duplicated(names(gene_list_plot))]

# 2. Gunakan fungsi secara spesifik dari library enrichplot
# Kita gunakan 'layout = "circular"' (sebagai string) bukan 'circular = TRUE' 
# untuk menghindari bentrokan dengan fungsi base R.
p_cnet <- enrichplot::cnetplot(ego_up, 
                               showCategory = 3, 
                               foldChange   = gene_list_plot, 
                               layout       = "circular") # Layout melingkar via ggraph

# 3. Tambahkan warna garis secara manual lewat ggplot2 (lebih aman)
p_cnet_final <- p_cnet + 
  ggplot2::labs(title = "Network Interaksi Gen-Jalur SLE",
                subtitle = "Gen melingkar untuk menghindari tumpukan teks",
                color = "Log2 FC") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.title = ggplot2::element_text(face="bold", size=14))

# 4. Tampilkan Plot
print(p_cnet)
ggsave("9_Cnetplot_Interaction_SLE.png", 
       plot = p_cnet_final, 
       width = 12, 
       height = 9, 
       dpi = 300)

cat("Selesai! File '9_Cnetplot_Interaction_SLE.png' telah siap untuk laporan Anda.\n")
