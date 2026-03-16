# analyze_BCR.R
# Purpose: Read BCR and TCR clonotype data from 10x VDJ output, summarize
#          clonotype sharing across timepoints, and visualize clonal overlap
#          for SLE patients treated with A-319.
# ==============================================================================

# ---- Load libraries ----------------------------------------------------------
library(tidyverse)

# ---- Set working directory ---------------------------------------------------
setwd("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/09_BCR_TCR")

# ==============================================================================
# SECTION 1: BCR Clonotype Analysis
# ==============================================================================

# ---- Read BCR clonotype files ------------------------------------------------
# Gather all VDJ-B clonotype CSVs across patients and timepoints.
b.clonotype.files <- list.files("TCR and BCR", recursive = TRUE,
                                all.files = TRUE, full.names = TRUE) %>%
  grep(pattern = "vdj_b", value = TRUE) %>%
  grep(pattern = "clonotype", value = TRUE)

b.clonotype.list <- list()
for (file in b.clonotype.files) {
  b.clonotype.list[[file]] <- read_csv(file)
  b.clonotype.list[[file]]$file.name <- file
}

b.clonotypes <- do.call(rbind, b.clonotype.list)

# ---- Parse patient and timepoint from file paths -----------------------------
b.clonotypes <- b.clonotypes %>%
  mutate(
    sample.timepoint = str_split_fixed(file.name, pattern = "/", n = 3)[, 2],
    patient          = str_split_fixed(sample.timepoint, pattern = " ", n = 2)[, 1],
    timepoint        = str_split_fixed(sample.timepoint, pattern = " ", n = 2)[, 2],
    clonotype_num    = gsub(x = clonotype_id, pattern = "clonotype",
                            replacement = "") %>% as.numeric()
  )

# ---- Visualize BCR clonotype counts and frequencies by sample ----------------
ggplot(b.clonotypes, aes(x = sample, fill = timepoint)) +
  geom_bar() + theme_bw()
ggplot(b.clonotypes, aes(x = sample, fill = timepoint, y = frequency)) +
  geom_col() + theme_bw()

# ---- Top clonotypes: proportion by sample and timepoint ----------------------
ggplot(b.clonotypes %>% filter(clonotype_num < 10),
       aes(x = clonotype_num, y = proportion)) +
  geom_point() + facet_grid(timepoint ~ sample, scales = "free")

# ---- Summarize BCR clonotype sharing across timepoints -----------------------
# For each patient, count how many timepoints each unique CDR3 appears in.
b.sum <- b.clonotypes %>%
  group_by(patient, cdr3s_aa) %>%
  summarize(
    frequency_sum          = sum(frequency, na.rm = TRUE),
    mean_prop              = mean(proportion),
    timepoints.represented = length(unique(timepoint)),
    base     = sum(timepoint == "Base"),
    D28      = sum(timepoint == "D28"),
    M3       = sum(timepoint == "M3"),
    timepoint = paste(unique(timepoint), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(timepoints.represented))

# ---- Compute per-patient clonotype overlap statistics (BCR) ------------------
n.clono.b <- b.sum %>%
  group_by(timepoints.represented, patient) %>%
  summarize(n.clono  = n(),
            max.prop = max(mean_prop),
            mean.prop = mean(mean_prop),
            n.cells  = sum(frequency_sum),
            .groups  = "drop") %>%
  left_join(
    b.sum %>%
      group_by(patient) %>%
      summarize(total.cells = sum(frequency_sum),
                total.clono = n(),
                .groups = "drop")
  ) %>%
  mutate(prop = n.cells / total.cells)

# ---- Summary: mean proportion of shared BCR clonotypes -----------------------
n.clono.b %>%
  group_by(timepoints.represented) %>%
  summarize(mean.prop = mean(prop), n = n(), sd.prop = sd(prop))

# ---- Compute per-patient clonotype overlap by timepoint combination ----------
n.clono.b2 <- b.sum %>%
  group_by(timepoint, patient) %>%
  summarize(n.clono  = n(),
            max.prop = max(mean_prop),
            mean.prop = mean(mean_prop),
            n.cells  = sum(frequency_sum),
            .groups  = "drop") %>%
  left_join(
    b.sum %>%
      group_by(patient) %>%
      summarize(total.cells = sum(frequency_sum),
                total.clono = n(),
                .groups = "drop")
  ) %>%
  mutate(prop = n.cells / total.cells)

# ---- Generate BCR overlap PDF ------------------------------------------------
pdf("BCR.overlap.pdf")

ggplot(b.sum %>% filter(timepoints.represented > 0),
       aes(x = timepoints.represented, y = frequency_sum, fill = timepoint)) +
  geom_col() + facet_wrap(~patient) + theme_bw() +
  ggtitle("BCR overlap") +
  scale_x_continuous(breaks = 1:3, limits = c(0.55, 3.5))

ggplot(n.clono.b,
       aes(x = patient, y = prop, fill = timepoints.represented, label = n.clono)) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5), color = "white") +
  geom_text(aes(label = total.cells, y = 1.02)) +
  theme_bw() +
  ggtitle("top: total cells, bars - middle: total clonotypes")

ggplot(n.clono.b2,
       aes(x = patient, y = prop, fill = timepoint, label = round(prop, 2))) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5), color = "white") +
  geom_text(aes(label = total.cells, y = 1.02)) +
  theme_bw() +
  ggtitle("top: total cells, bars - middle: proportion")

dev.off()

# ==============================================================================
# SECTION 2: TCR Clonotype Analysis
# ==============================================================================

# ---- Read TCR clonotype files ------------------------------------------------
t.clonotype.files <- list.files("TCR and BCR", recursive = TRUE,
                                all.files = TRUE, full.names = TRUE) %>%
  grep(pattern = "vdj_t", value = TRUE) %>%
  grep(pattern = "clonotype", value = TRUE)

t.clonotype.list <- list()
for (file in t.clonotype.files) {
  t.clonotype.list[[file]] <- read_csv(file)
  t.clonotype.list[[file]]$file.name <- file
}

t.clonotypes <- do.call(rbind, t.clonotype.list)

# ---- Parse patient and timepoint from file paths -----------------------------
t.clonotypes <- t.clonotypes %>%
  mutate(
    sample.timepoint = str_split_fixed(file.name, pattern = "/", n = 3)[, 2],
    patient          = str_split_fixed(sample.timepoint, pattern = " ", n = 2)[, 1],
    sample           = patient,
    timepoint        = str_split_fixed(sample.timepoint, pattern = " ", n = 2)[, 2],
    clonotype_num    = gsub(x = clonotype_id, pattern = "clonotype",
                            replacement = "") %>% as.numeric()
  )

# ---- Visualize TCR clonotype counts and frequencies by sample ----------------
ggplot(t.clonotypes, aes(x = sample, fill = timepoint)) +
  geom_bar() + theme_bw()
ggplot(t.clonotypes, aes(x = sample, fill = timepoint, y = frequency)) +
  geom_col() + theme_bw()
ggplot(t.clonotypes %>% filter(clonotype_num < 50),
       aes(x = clonotype_num, y = frequency)) +
  geom_point() + facet_grid(timepoint ~ sample, scales = "free")

# ---- Summarize TCR clonotype sharing across timepoints -----------------------
t.sum <- t.clonotypes %>%
  group_by(patient, cdr3s_aa) %>%
  summarize(
    frequency_sum          = sum(frequency, na.rm = TRUE),
    mean_prop              = mean(proportion),
    timepoints.represented = length(unique(timepoint)),
    base     = sum(timepoint == "Base"),
    D28      = sum(timepoint == "D28"),
    M3       = sum(timepoint == "M3"),
    timepoint = paste(unique(timepoint), collapse = ","),
    .groups = "drop"
  ) %>%
  arrange(desc(timepoints.represented), desc(frequency_sum))

t.sum$cdr3s_aa <- factor(t.sum$cdr3s_aa, levels = unique(t.sum$cdr3s_aa))

# ---- Compute per-patient clonotype overlap statistics (TCR) ------------------
n.clono.t <- t.sum %>%
  group_by(timepoints.represented, patient) %>%
  summarize(n.clono  = n(),
            max.prop = max(mean_prop),
            mean.prop = mean(mean_prop),
            n.cells  = sum(frequency_sum),
            .groups  = "drop") %>%
  left_join(
    t.sum %>%
      group_by(patient) %>%
      summarize(total.cells = sum(frequency_sum),
                total.clono = n(),
                .groups = "drop")
  ) %>%
  mutate(prop = n.cells / total.cells)

n.clono.t %>%
  group_by(timepoints.represented) %>%
  summarize(mean.prop = mean(prop), n = n(), sd.prop = sd(prop))

# ---- Generate TCR overlap PDF ------------------------------------------------
pdf("TCR.overlap.pdf")

ggplot(t.sum %>% filter(timepoints.represented > 0),
       aes(x = timepoints.represented, y = frequency_sum)) +
  geom_col(aes(fill = timepoint)) + facet_wrap(~patient) +
  theme_bw() + ggtitle("TCR overlap")

ggplot(t.sum %>% filter(timepoints.represented > 0),
       aes(x = timepoints.represented, y = frequency_sum)) +
  geom_col(aes(fill = timepoint)) + facet_wrap(~patient) +
  theme_bw() + ggtitle("TCR overlap, (mean prop)") +
  geom_text(aes(x = timepoints.represented, y = 5000,
                label = paste0(n.clono, "\\n (", round(mean.prop * 100, 1), "%)")),
            data = n.clono.t)

ggplot(t.sum %>% filter(timepoints.represented > 0),
       aes(x = timepoints.represented, y = frequency_sum)) +
  geom_col(aes(fill = timepoint)) + facet_wrap(~patient) +
  theme_bw() + ggtitle("TCR overlap, (max prop)") +
  geom_text(aes(x = timepoints.represented, y = 5000,
                label = paste0(n.clono, "\\n (", round(max.prop * 100, 1), "%)")),
            data = n.clono.t)

ggplot(n.clono.t,
       aes(x = patient, y = prop, fill = timepoints.represented, label = n.clono)) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5), color = "white") +
  theme_bw() +
  ggtitle("top: total cells, bars - middle: total clonotypes")

ggplot(n.clono.t,
       aes(x = patient, y = prop, fill = timepoints.represented, label = round(prop, 2))) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5), color = "white") +
  geom_text(aes(label = total.cells, y = 1.02)) +
  theme_bw() +
  ggtitle("top: total cells, bars - middle: proportion")

dev.off()

# ==============================================================================
# SECTION 3: Aggregated TCR / BCR clonotype count summaries
# ==============================================================================

pdf("t.b.clonotype.count.pdf")

# ---- Grouped TCR summary: by number of timepoints represented ---------------
n.clono.t.group <- t.sum %>%
  group_by(timepoints.represented) %>%
  summarize(n.clono  = n(),
            max.prop = max(mean_prop),
            mean.prop = mean(mean_prop),
            n.cells  = sum(frequency_sum),
            .groups  = "drop") %>%
  mutate(prop = n.cells / 96292)

ggplot(n.clono.t.group,
       aes(x = "", y = prop, fill = timepoints.represented,
           label = paste0(n.clono, " (n=", n.cells, ")"))) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5), color = "white") +
  theme_bw() + ggtitle("TCR: clonotype overlap by # shared timepoints")

ggplot(n.clono.t.group,
       aes(x = "", y = n.cells, fill = timepoints.represented,
           label = paste0(n.clono, " (n=", n.cells, ")"))) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5), color = "white") +
  theme_bw() + ggtitle("TCR: cell counts by # shared timepoints")

# ---- Grouped TCR summary: by individual timepoint ----------------------------
n.clono.t.group2 <- t.sum %>%
  group_by(timepoint) %>%
  summarize(n.clono  = n(),
            max.prop = max(mean_prop),
            mean.prop = mean(mean_prop),
            n.cells  = sum(frequency_sum),
            .groups  = "drop") %>%
  mutate(prop = n.cells / 96292)

ggplot(n.clono.t.group2,
       aes(x = "", y = prop, fill = timepoint,
           label = paste0(n.clono, " (n=", n.cells, ")"))) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5), color = "white") +
  theme_bw() + ggtitle("TCR: clonotype proportions by timepoint")

ggplot(n.clono.t.group2,
       aes(x = "", y = n.cells, fill = timepoint,
           label = paste0(n.clono, " (n=", n.cells, ")"))) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5), color = "white") +
  theme_bw() + ggtitle("TCR: cell counts by timepoint")

# ---- Grouped BCR summary: by number of timepoints represented ---------------
n.clono.b.group <- b.sum %>%
  group_by(timepoints.represented) %>%
  summarize(n.clono  = n(),
            max.prop = max(mean_prop),
            mean.prop = mean(mean_prop),
            n.cells  = sum(frequency_sum),
            .groups  = "drop") %>%
  mutate(prop = n.cells / 5953)

ggplot(n.clono.b.group,
       aes(x = "", y = prop, fill = timepoints.represented,
           label = paste0(n.clono, " (n=", n.cells, ")"))) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5), color = "white") +
  theme_bw() + ggtitle("BCR: clonotype overlap by # shared timepoints")

ggplot(n.clono.b.group,
       aes(x = "", y = n.cells, fill = timepoints.represented,
           label = paste0(n.clono, " (n=", n.cells, ")"))) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5), color = "white") +
  theme_bw() + ggtitle("BCR: cell counts by # shared timepoints")

# ---- Grouped BCR summary: by individual timepoint ----------------------------
n.clono.b.group2 <- b.sum %>%
  group_by(timepoint) %>%
  summarize(n.clono  = n(),
            max.prop = max(mean_prop),
            mean.prop = mean(mean_prop),
            n.cells  = sum(frequency_sum),
            .groups  = "drop") %>%
  mutate(prop = n.cells / 5953)

ggplot(n.clono.b.group2,
       aes(x = "", y = prop, fill = timepoint,
           label = paste0(n.clono, " (n=", n.cells, ")"))) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5), color = "white") +
  theme_bw() + ggtitle("BCR: clonotype proportions by timepoint")

ggplot(n.clono.b.group2,
       aes(x = "", y = n.cells, fill = timepoint,
           label = paste0(n.clono, " (n=", n.cells, ")"))) +
  geom_col() +
  geom_text(position = position_stack(vjust = 0.5), color = "white") +
  theme_bw() + ggtitle("BCR: cell counts by timepoint")

dev.off()

# ==============================================================================
# SECTION 4: Parse chain-level information for TCR and BCR
# ==============================================================================

# ---- TCR chain parsing -------------------------------------------------------
t.clonotypes <- t.clonotypes %>%
  mutate(
    chain.1       = str_split_fixed(cdr3s_aa, pattern = ";", n = 2)[, 1],
    chain.2       = str_split_fixed(cdr3s_aa, pattern = ";", n = 2)[, 2],
    chain.1.short = str_split_fixed(cdr3s_aa, pattern = ":", n = 2)[, 1],
    chain.2.short = str_split_fixed(chain.2,  pattern = ":", n = 2)[, 1]
  )

# ---- BCR chain parsing -------------------------------------------------------
b.clonotypes <- b.clonotypes %>%
  mutate(
    chain.1       = str_split_fixed(cdr3s_aa, pattern = ";", n = 2)[, 1],
    chain.2       = str_split_fixed(cdr3s_aa, pattern = ";", n = 2)[, 2],
    chain.1.short = str_split_fixed(cdr3s_aa, pattern = ":", n = 2)[, 1],
    chain.2.short = str_split_fixed(chain.2,  pattern = ":", n = 2)[, 1]
  )

# ---- Save results ------------------------------------------------------------
save.image("2025_10_20_image_tb_clonotypes.rdata")
saveRDS(t.clonotypes, "t.clonotypes.rds")
saveRDS(b.clonotypes, "b.clonotypes.rds")