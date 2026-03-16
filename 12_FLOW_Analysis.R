setwd("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/28_B_cell_subsets_flow")

B_cell_subsets_tidy <- readRDS("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/28_B_cell_subsets_flow/B_cell_subsets_tidy.rds")
B_cell_subsets_tidy
bcell_ratios <- readRDS("~/Library/CloudStorage/Box-Box/Tan_Lab/Collaborations/XQ_Yan/scRNA/28_B_cell_subsets_flow/all_files_20260207T223900/bcell_ratios.rds")

ggplot(B_cell_subsets_tidy, aes(x = day, y = value)) + geom_point() + geom_line(aes(group = patient)) + facet_wrap(~b_cell_subset) + geom_vline(xintercept = 28, linetype = 'dashed', color = "red")


ggplot(B_cell_subsets_tidy %>% filter(day %in% c(1, 3, 8, 10, 15, 22, 26, 28, 42, 60, 90, 120, 150, 180, 270, 360)), aes(x = day, y = value)) + geom_point() + geom_line(aes(group = patient)) + facet_wrap(~b_cell_subset) + geom_vline(xintercept = 28, linetype = 'dashed', color = "red")


ggplot(B_cell_subsets_tidy %>% filter(day %in% c(1, 26, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = day, y = value)) + geom_point() + geom_line(aes(group = patient)) + facet_wrap(~b_cell_subset) + geom_vline(xintercept = 28, linetype = 'dashed', color = "red") + theme_bw()



ggplot(B_cell_subsets_tidy %>% filter(day %in% c(1, 26, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = day, y = value, color = dose)) + geom_point() + geom_line(aes(group = patient)) + facet_wrap(~b_cell_subset) + 
  geom_vline(xintercept = 28, linetype = 'dashed', color = "red") + theme_bw() + scale_x_continuous(breaks = c(1, 26, 42, 60, 90, 120, 150, 180, 270, 360))



B_cell_subsets_tidy$timepoint_label_factor = factor(B_cell_subsets_tidy$timepoint_label, levels = B_cell_subsets_tidy$timepoint_label %>% unique())

ggplot(B_cell_subsets_tidy %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = timepoint_label_factor, y = value, color = dose)) + geom_point() + geom_line(aes(group = patient)) + facet_wrap(~b_cell_subset) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw()

ggplot(B_cell_subsets_tidy %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = timepoint_label_factor, y = fold_change_from_baseline, color = dose)) + geom_point() + geom_line(aes(group = patient)) + facet_wrap(~b_cell_subset) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw() + scale_y_log10()



B_cell_subsets_tidy$b_cell_subset = factor(B_cell_subsets_tidy$b_cell_subset, levels = c("naive B", "Switched memory B", "Unswitched memory B", "Plasmablast"))

dose.summary.b.cells = B_cell_subsets_tidy %>% group_by(timepoint_label_factor, b_cell_subset, dose, day) %>% summarize(mean.fc = mean(fold_change_from_baseline, na.rm = T), 
                                                                                            max.fc = max(fold_change_from_baseline, na.rm = T), 
                                                                                            min.fc = min(fold_change_from_baseline, na.rm = T),
                                                                                            mean.value = mean(value, na.rm = T), 
                                                                                            sd = sd(value, na.rm = T), 
                                                                                            max.val = max(value, na.rm = T), 
                                                                                            min.val = min(value, na.rm = T))

dose.summary.b.cells = dose.summary.b.cells %>% filter(b_cell_subset != "IgD-CD27- DN B")
B_cell_subsets_tidy = B_cell_subsets_tidy %>% filter(b_cell_subset != "IgD-CD27- DN B")


pdf("B.cell.overview.pdf")
ggplot(dose.summary.b.cells %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = timepoint_label_factor, y = mean.value, color = dose)) + geom_point() + geom_line(aes(group = dose), linewidth = 2) + facet_wrap(~b_cell_subset) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw() +   
  geom_line(data = B_cell_subsets_tidy %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
                                                  aes(group = patient, x = timepoint_label_factor, y = value, color = dose), alpha = 0.5)


ggplot(dose.summary.b.cells %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = timepoint_label_factor, y = mean.value, color = dose)) + geom_point() + geom_line(aes(group = dose), linewidth = 2) + facet_wrap(~b_cell_subset) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw() +   
  geom_errorbar(aes(ymin = min.val, ymax = max.val), width = 0, size = 0.5, alpha = 0.7, position =position_jitter(width = 0.2)) + ggtitle("min.max")


ggplot(dose.summary.b.cells %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = timepoint_label_factor, y = mean.value, color = dose)) + geom_point() + geom_line(aes(group = dose), linewidth = 2) + facet_wrap(~b_cell_subset) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw() +   
  geom_errorbar(aes(ymin = mean.value - sd, ymax = mean.value + sd), width = 0, size = 0.5, alpha = 0.7, position =position_jitter(width = 0.1)) + ggtitle("sd")



ggplot(dose.summary.b.cells %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = timepoint_label_factor, y = mean.fc, color = dose)) + geom_point() + geom_line(aes(group = dose), linewidth = 2) + facet_wrap(~b_cell_subset) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw() +   
  geom_errorbar(aes(ymin = min.fc, ymax = max.fc), width = 0, size = 0.5, alpha = 0.7, position =position_jitter(width = 0.2)) + ggtitle("min.max FC") + scale_y_log10()

ggplot(dose.summary.b.cells.all %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = timepoint_label_factor, y = mean.value, color = b_cell_subset)) + geom_point() + geom_line(aes(group = b_cell_subset), linewidth = 2) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw() +
  geom_errorbar(aes(ymin = min.val, ymax = max.val), width = 0, size = 0.5, alpha = 0.7, position =position_jitter(width = 0.03, seed = 42)) + ggtitle("min.max")

ggplot(dose.summary.b.cells.all %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360), b_cell_subset != "IgD-CD27- DN B"), 
       aes(x = timepoint_label_factor, y = mean.value, color = b_cell_subset)) + geom_point() + geom_line(aes(group = b_cell_subset), linewidth = 2) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw() +
  geom_errorbar(aes(ymin = min.val, ymax = max.val), width = 0, size = 0.5, alpha = 0.7, position =position_jitter(width = 0.03, seed = 42)) + ggtitle("min.max")


dev.off()


dose.summary.b.cells.all = B_cell_subsets_tidy %>% group_by(timepoint_label_factor, b_cell_subset, day) %>% summarize(mean.fc = mean(fold_change_from_baseline, na.rm = T), 
                                                                                                                        max.fc = max(fold_change_from_baseline, na.rm = T), 
                                                                                                                        min.fc = min(fold_change_from_baseline, na.rm = T),
                                                                                                                        mean.value = mean(value, na.rm = T), 
                                                                                                                        sd = sd(value, na.rm = T), 
                                                                                                                        max.val = max(value, na.rm = T), 
                                                                                                                        min.val = min(value, na.rm = T))

pdf("Naive_memory.pdf", height = 6, width = 8)
ggplot(dose.summary.b.cells %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360), b_cell_subset %in% c("naive B", "Switched memory B")), 
       aes(x = timepoint_label_factor, y = mean.value, color = b_cell_subset)) + geom_point() + geom_line(aes(group = b_cell_subset), linewidth = 2) + facet_wrap(~dose) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw()

ggplot(dose.summary.b.cells.all %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360), b_cell_subset %in% c("naive B", "Switched memory B")), 
       aes(x = timepoint_label_factor, y = mean.value, color = b_cell_subset)) + geom_point() + geom_line(aes(group = b_cell_subset), linewidth = 2) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw() +
  geom_errorbar(aes(ymin = min.val, ymax = max.val), width = 0, size = 0.5, alpha = 0.7, position =position_jitter(width = 0.03, seed = 42)) + ggtitle("min.max")


ggplot(dose.summary.b.cells.all %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360), b_cell_subset %in% c("naive B", "Switched memory B")), 
       aes(x = timepoint_label_factor, y = mean.value, color = b_cell_subset)) + geom_point() + geom_line(aes(group = b_cell_subset), linewidth = 2) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw() +
  geom_errorbar(aes(ymin = mean.value-sd, ymax = mean.value+sd), width = 0, size = 0.5, alpha = 0.7, position =position_jitter(width = 0.03, seed = 42)) + ggtitle("sd")

ggplot(dose.summary.b.cells.all %>% filter(day %in% c(1, 8, 15, 22, 28, 42, 60, 90, 120, 150, 180, 270, 360), b_cell_subset %in% c("naive B", "Switched memory B")), 
       aes(x = timepoint_label_factor, y = mean.value, color = b_cell_subset)) + geom_point() + geom_line(aes(group = b_cell_subset), linewidth = 3) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw() +
  geom_line(data = B_cell_subsets_tidy %>% filter(day %in% c(1, 8, 15, 22, 28, 42, 60, 90, 120, 150, 180, 270, 360),  b_cell_subset %in% c("Switched memory B")), 
            aes(group = patient, x = timepoint_label_factor, y = value, color = b_cell_subset), alpha = 0.5, linewidth = 0.5)+
  geom_line(data = B_cell_subsets_tidy %>% filter(day %in% c(1, 8, 15, 22, 28, 42, 60, 90, 120, 150, 180, 270, 360),  b_cell_subset %in% c("naive B")), 
            aes(group = patient, x = timepoint_label_factor, y = value, color = b_cell_subset), alpha = 0.5, linewidth = 0.5)


dev.off()







dose.summary.b.cells2 = bcell_ratios %>% group_by(Analysis_Timepoint, Dose_Group, Day_Numeric) %>% summarize(mean.ratio = mean(Naive_Memory_Ratio, na.rm = T), 
                                                                                                                        max.ratio = max(Naive_Memory_Ratio, na.rm = T), 
                                                                                                                        min.ratio = min(Naive_Memory_Ratio, na.rm = T),
                                                                                                                        sd = sd(Naive_Memory_Ratio, na.rm = T))


pdf("NBC_MBC_ratio.pdf", height = 6, width = 8)
ggplot(bcell_ratios %>% filter(Day_Numeric %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = Analysis_Timepoint, y = Log2_Naive_Memory_Ratio, color = Dose_Group)) + geom_point() + geom_line(aes(group = Patient)) +
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw()



ggplot(bcell_ratios %>% filter(Day_Numeric %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = Analysis_Timepoint, y = log(Naive_Memory_Ratio, base = 10), color = Dose_Group)) + geom_point() + geom_line(aes(group = Patient)) +
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw()

ggplot(dose.summary.b.cells2 %>% filter(Day_Numeric %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = Analysis_Timepoint, y = log(mean.ratio, base = 10), color = Dose_Group)) + geom_point() + geom_line(aes(group = Dose_Group)) +
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw() +   
  geom_errorbar(aes(ymin = log(min.ratio, base = 10), ymax = log(max.ratio, base = 10)), width = 0, size = 0.5, alpha = 0.7, position =position_jitter(width = 0.2)) + ggtitle("Ratio of NBC / MBC, min.max")



ggplot(bcell_ratios %>% filter(Day_Numeric %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = Analysis_Timepoint, y = log(Naive_Memory_Ratio, base = 10))) + geom_boxplot(outlier.size = -1) +  
  #geom_point(aes(color = Dose_Group, shape = Dose_Group), position = position_jitter(width = 0.1))  + 
  theme_bw()+
  geom_line(data = dose.summary.b.cells2%>% filter(Day_Numeric %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
            aes(x = Analysis_Timepoint, y = log(mean.ratio, base = 10), color = Dose_Group, group = Dose_Group),linewidth = 2) + stat_compare_means() +
  geom_errorbar(data = dose.summary.b.cells2%>% filter(Day_Numeric %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
                aes(x = Analysis_Timepoint, y = log(mean.ratio, base = 10), color = Dose_Group, group = Dose_Group, ymin = log(min.ratio, base = 10), 
                    ymax = log(max.ratio, base = 10)), width = 0, size = 0.5, alpha = 0.7, position =position_jitter(width = 0.05)) + ggtitle("Ratio of NBC / MBC, min.max")



ggplot(bcell_ratios %>% filter(Day_Numeric %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = Analysis_Timepoint, y = log(Naive_Memory_Ratio, base = 10))) + 
  #geom_point(aes(color = Dose_Group, shape = Dose_Group), position = position_jitter(width = 0.1))  + 
  theme_bw()+
  geom_line(data = dose.summary.b.cells2%>% filter(Day_Numeric %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
            aes(x = Analysis_Timepoint, y = log(mean.ratio, base = 10), color = Dose_Group, group = Dose_Group),linewidth = 2) + stat_compare_means() +
  geom_boxplot(outlier.size = -1, alpha = 0.5, fill = NA) 

ggplot(bcell_ratios %>% filter(Day_Numeric %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = Analysis_Timepoint, y = log(Naive_Memory_Ratio, base = 10))) + 
  #geom_boxplot(outlier.size = -1) +  
  geom_line(aes(group = Patient, color = Dose_Group), alpha = 0.3) +
  geom_point(aes(color = Dose_Group, shape = Dose_Group))  + 
  theme_bw()+
  geom_line(data = dose.summary.b.cells2%>% filter(Day_Numeric %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
            aes(x = Analysis_Timepoint, y = log(mean.ratio, base = 10), color = Dose_Group, group = Dose_Group),linewidth = 2) + stat_compare_means()
dev.off()



#patient 5

B_cell_subsets_tidy = B_cell_subsets_tidy %>% mutate(p5.comp = case_when(patient %in% c("P1", "P2", "P4", "P3", "P6")~ "Other C1",
                                                    patient == "P5" ~ "P5",
                                                   patient %in% c("P6", "P7", "P8", "P9", "P10", "P11", "P12") ~ "C2.C3"))

dose.summary.b.cells.p5 = B_cell_subsets_tidy %>% group_by(timepoint_label_factor, b_cell_subset, day, patient == "P5") %>% summarize(mean.fc = mean(fold_change_from_baseline, na.rm = T), 
                                                                                                                      max.fc = max(fold_change_from_baseline, na.rm = T), 
                                                                                                                      min.fc = min(fold_change_from_baseline, na.rm = T),
                                                                                                                      mean.value = mean(value, na.rm = T), 
                                                                                                                      sd = sd(value, na.rm = T), 
                                                                                                                      max.val = max(value, na.rm = T), 
                                                                                                                      min.val = min(value, na.rm = T))

dose.summary.b.cells.p5.2 = B_cell_subsets_tidy %>% group_by(timepoint_label_factor, b_cell_subset, day, p5.comp) %>% summarize(mean.fc = mean(fold_change_from_baseline, na.rm = T), 
                                                                                                                                      max.fc = max(fold_change_from_baseline, na.rm = T), 
                                                                                                                                      min.fc = min(fold_change_from_baseline, na.rm = T),
                                                                                                                                      mean.value = mean(value, na.rm = T), 
                                                                                                                                      sd = sd(value, na.rm = T), 
                                                                                                                                      max.val = max(value, na.rm = T), 
                                                                                                                                      min.val = min(value, na.rm = T))


bcell_ratios = bcell_ratios %>% mutate(p5.comp = case_when(Patient %in% c("P1", "P2", "P4", "P3", "P6")~ "Other C1",
                                                            Patient == "P5" ~ "P5",
                                                             Patient %in% c("P6", "P7", "P8", "P9", "P10", "P11", "P12") ~ "C2.C3"))

bcell_ratios.p5 = bcell_ratios %>% group_by(Analysis_Timepoint, Patient == "P5", Day_Numeric) %>% summarize(mean.ratio = mean(Naive_Memory_Ratio, na.rm = T), 
                                                                                                             max.ratio = max(Naive_Memory_Ratio, na.rm = T), 
                                                                                                             min.ratio = min(Naive_Memory_Ratio, na.rm = T),
                                                                                                             sd = sd(Naive_Memory_Ratio, na.rm = T))

bcell_ratios.p5.2 = bcell_ratios %>% group_by(Analysis_Timepoint, p5.comp, Day_Numeric) %>% summarize(mean.ratio = mean(Naive_Memory_Ratio, na.rm = T), 
                                                                                                            max.ratio = max(Naive_Memory_Ratio, na.rm = T), 
                                                                                                            min.ratio = min(Naive_Memory_Ratio, na.rm = T),
                                                                                                            sd = sd(Naive_Memory_Ratio, na.rm = T))

pdf("p5.outlier.pdf", height = 6, width = 8)
ggplot(dose.summary.b.cells.p5 %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = timepoint_label_factor, y = mean.value, color =  `patient == "P5"`)) + geom_point() + geom_line(aes(group = `patient == "P5"`), linewidth = 2) + facet_wrap(~b_cell_subset) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw()

ggplot(dose.summary.b.cells.p5.2 %>% filter(day %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = timepoint_label_factor, y = mean.value, color =  p5.comp)) + geom_point() + geom_line(aes(group = p5.comp), linewidth = 2) + facet_wrap(~b_cell_subset) + 
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw()


ggplot(bcell_ratios.p5 %>% filter(Day_Numeric %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = Analysis_Timepoint, y = log(mean.ratio, base = 10), color = `Patient == "P5"`)) + geom_point() + geom_line(aes(group = `Patient == "P5"`)) +
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw()

ggplot(bcell_ratios.p5.2 %>% filter(Day_Numeric %in% c(1, 28, 42, 60, 90, 120, 150, 180, 270, 360)), 
       aes(x = Analysis_Timepoint, y = log(mean.ratio, base = 10), color = p5.comp)) + geom_point() + geom_line(aes(group = p5.comp)) +
  geom_vline(xintercept = "D28", linetype = 'dashed', color = "red") + theme_bw()


dev.off()


save.image("B_cell_subsets_CRS_25_28.rdata")



days.to.loop =  c(1, 26, 42, 60, 90, 120, 150, 180, 270, 360)

for(x in days.to.loop){
  print(x)
  to.test1 = B_cell_subsets_tidy %>% filter(day == x, b_cell_subset %in% c("Switched memory B"), !is.na(value)) %>% as_tibble()
  to.test2 = B_cell_subsets_tidy %>% filter(day == x, b_cell_subset %in% c("naive B"), !is.na(value)) %>% as_tibble()
  t.test(to.test1$value, to.test2$value, paired = T) %>% print()
  
  
}
