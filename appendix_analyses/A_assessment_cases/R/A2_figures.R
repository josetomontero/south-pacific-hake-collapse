# Figures of the assessment-case sensitivity (A1_fit_P3_cases.R).
# A2_cases_biomass: total biomass 1997-2021 of Cases 0-3 and the series of the paper.
# A2_cases_dAICc: dAICc of the 8 candidate models of 2001-2020 fitted to each series (U = official
#   catch / case biomass); the squid predation model (8) highlighted.
A <- rprojroot::find_root(rprojroot::has_file("R/A1_fit_P3_cases.R"))
PROJ <- A; RES <- A
source(file.path(A, "R", "lib", "fig_theme.R"))
HAKE <- normalizePath(file.path(A, "..", ".."))   # repository root 

hd <- read.delim(file.path(HAKE, "data", "derived", "hake_dataset.txt"))
cs <- read.csv(file.path(A, "data", "IFOP2022_cases0-3_1997_2021.csv"))
tab <- read.delim(file.path(A, "tables", "A1_P3_models_by_case.txt"))

SER_COLS <- c("Paper series" = "black", "Case 0" = "#D55E00", "Case 1" = "#B8860B", "Case 2" = "#0072B2", "Case 3" = "#009E73")
b <- rbind(data.frame(series = paste("Case", cs$case), year = cs$year, B = cs$total_biomass_kt),
           data.frame(series = "Paper series", year = hd$year, B = hd$biomass / 1000))
b <- subset(b, year >= 1997 & year <= 2021 & !is.na(B))
b$series <- factor(b$series, names(SER_COLS))
p1 <- ggplot(b, aes(year, B, colour = series, linewidth = series == "Paper series")) +
  geom_line() + geom_point(size = 0.8) +
  scale_colour_manual(values = SER_COLS,
                      labels = c("Paper series", "Case 0: official landings", "Case 1: official landings, Francis weighting",
                                 "Case 2: CCT catches", "Case 3: CCT catches, Francis weighting")) +
  scale_linewidth_manual(values = c(`TRUE` = 1, `FALSE` = 0.5), guide = "none") +
  scale_x_continuous(breaks = seq(1997, 2021, 4)) + scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.05))) +
  labs(x = "Year", y = "Total biomass 2+ (kt)") + theme_pub() + theme(legend.position = c(0.68, 0.8))
save_fig(p1, "A2_cases_biomass", width_cm = 12, height_cm = 8, dir = file.path(A, "figs"))

d <- subset(tab, !grepl("\\[u of the case table\\]", series))
d$series <- factor(sub(":.*", "", d$series), names(SER_COLS))
d$model <- factor(d$model, rev(unique(d$model)))
p2 <- ggplot(d, aes(dAICc, model, fill = model == "8 squid predation")) +
  geom_col(width = 0.7) + geom_vline(xintercept = 2, linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ series, nrow = 1) +
  scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "grey60"), guide = "none") +
  labs(x = expression(Delta * "AICc (2001-2020 models)"), y = NULL) + theme_pub() +
  theme(strip.background = element_blank(), strip.text = element_text(size = 8))
save_fig(p2, "A2_cases_dAICc", width_cm = 18, height_cm = 7, dir = file.path(A, "figs"))
