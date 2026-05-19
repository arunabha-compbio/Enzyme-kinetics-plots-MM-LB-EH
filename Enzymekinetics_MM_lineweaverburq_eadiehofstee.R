# Load required libraries
library(ggplot2)

# ----- Experimental Constants -----
epsilon         <- 18800      # Molar extinction coefficient of pNP at 405 nm (M⁻¹ cm⁻¹)
path_length     <- 1          # Cuvette path length in cm
reaction_time   <- 30         # *** Set your actual reaction time in minutes ***

# Beer-Lambert conversion factor: OD → mM/min
# V (mM/min) = OD × 1000 / (ε × l × t)
# Therefore: 1/V = (ε × l × t) / (OD × 1000)
BL_factor <- (epsilon * path_length * reaction_time) / 1000

cat("=== Beer-Lambert Conversion ===\n")
cat(sprintf("ε (pNP, 405nm) : %d M⁻¹cm⁻¹\n", epsilon))
cat(sprintf("Path length    : %d cm\n", path_length))
cat(sprintf("Reaction time  : %d min\n", reaction_time))
cat(sprintf("Scaling factor : %.2f  [1/V = factor × 1/OD]\n\n", BL_factor))

# ----- Raw Data -----
inv_S  <- c(1, 1.25, 1.67, 2.5, 5)         # 1/[S] in mM⁻¹
inv_OD <- c(2.7, 3.22, 3.85, 5.26, 9.09)   # 1/OD

# ----- Convert 1/OD → 1/V (min/mM) -----
# 1/V = (ε × l × t / 1000) × (1/OD)
inv_V  <- inv_OD * BL_factor

df <- data.frame(inv_S, inv_OD, inv_V)

cat("=== Converted Data ===\n")
print(data.frame(
  "1/[S] (mM⁻¹)" = inv_S,
  "1/OD"          = inv_OD,
  "1/V (min/mM)"  = round(inv_V, 2),
  check.names     = FALSE
))
cat("\n")

# ----- Linear Regression on converted values -----
lm_fit    <- lm(inv_V ~ inv_S, data = df)
intercept <- coef(lm_fit)[1]   # 1/Vmax
slope     <- coef(lm_fit)[2]   # Km/Vmax

# ----- Kinetic Parameters -----
Vmax  <- 1 / intercept         # mM/min
Km    <- slope * Vmax          # mM
x_int <- -intercept / slope    # -1/Km

cat("=== Lineweaver-Burk Kinetic Parameters ===\n")
cat(sprintf("Slope       : %.4f\n", slope))
cat(sprintf("Y-intercept : %.4f min/mM  →  Vmax = %.4f mM/min\n", intercept, Vmax))
cat(sprintf("X-intercept : %.4f mM⁻¹   →  Km   = %.4f mM\n", x_int, Km))
cat(sprintf("R²          : %.4f\n", summary(lm_fit)$r.squared))

# ----- Extended regression line -----
x_range  <- seq(x_int - 0.3, max(inv_S) + 0.2, length.out = 500)
fit_line <- data.frame(
  inv_S = x_range,
  inv_V = intercept + slope * x_range
)

# ----- Intercept points -----
y_int_pt <- data.frame(inv_S = 0,     inv_V = intercept)
x_int_pt <- data.frame(inv_S = x_int, inv_V = 0)

# ----- Plot -----
ggplot() +
  
  # Regression line
  geom_line(data = fit_line,
            aes(x = inv_S, y = inv_V),
            color = "#E05A2B", linewidth = 1.2) +
  
  # Reference axes
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.5) +
  
  # Data points
  geom_point(data = df,
             aes(x = inv_S, y = inv_V),
             color = "#2B4BA8", size = 4, shape = 21,
             fill = "#5B8DEF", stroke = 1.2) +
  
  # Y-intercept (1/Vmax)
  geom_point(data = y_int_pt,
             aes(x = inv_S, y = inv_V),
             color = "#228B22", size = 4, shape = 17) +
  annotate("text",
           x = 0.08, y = intercept - 8,
           label = paste0("1/Vmax = ", round(intercept, 2), " min/mM",
                          "\nVmax = ",  round(Vmax, 4), " mM/min"),
           hjust = 0, color = "#228B22", size = 3.5, fontface = "italic") +
  
  # X-intercept (-1/Km)
  geom_point(data = x_int_pt,
             aes(x = inv_S, y = inv_V),
             color = "#8B0000", size = 4, shape = 17) +
  annotate("text",
           x = x_int + 0.08, y = 15,
           label = paste0("-1/Km = ", round(x_int, 3), " mM⁻¹",
                          "\nKm = ",   round(Km, 4), " mM"),
           hjust = 0, color = "#8B0000", size = 3.5, fontface = "italic") +
  
  # R² annotation
  annotate("text",
           x = 3.0, y = 280,
           label = paste0("R² = ", round(summary(lm_fit)$r.squared, 4)),
           hjust = 0, color = "grey30", size = 3.8, fontface = "italic") +
  
  # Labels
  labs(
    title    = "Lineweaver–Burk Plot: Alkaline Phosphatase Kinetics",
    subtitle = paste0("pNPP substrate  |  ε = 18,800 M⁻¹cm⁻¹  |  t = ",
                      reaction_time, " min  |  Axes in true kinetic units"),
    x        = "1/[S]  (mM⁻¹)",
    y        = "1/V  (min / mM)"
  ) +
  
  scale_x_continuous(breaks = seq(-1, 6, by = 1)) +
  
  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey50", size = 9.5),
    axis.title       = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey92")
  )
# ============================================================
# ENZYME KINETICS ANALYSIS — ALKALINE PHOSPHATASE (pNPP)
# Plots: Michaelis-Menten & Eadie-Hofstee
# Kinetic parameters sourced from Lineweaver-Burk regression
# ============================================================

library(ggplot2)

# ----- Experimental Constants -----
epsilon       <- 18800   # Molar extinction coefficient of pNP at 405nm (M⁻¹cm⁻¹)
path_length   <- 1       # cm
reaction_time <- 30      # minutes — change to your actual assay time

# Beer-Lambert conversion: V (mM/min) = OD × 1000 / (ε × l × t)
BL_factor <- (epsilon * path_length * reaction_time) / 1000

# ----- Raw Data -----
S_mM     <- c(1.0, 0.8, 0.6, 0.4, 0.2)      # [S] in mM
od_raw   <- c(0.37, 0.31, 0.26, 0.19, 0.11)  # OD at 405nm

# ----- Convert OD → V (mM/min) via Beer-Lambert -----
V_mM_min <- od_raw / BL_factor

df <- data.frame(S = S_mM, V = V_mM_min)

cat("=== Converted Velocity Values ===\n")
print(data.frame(
  "[S] (mM)"   = S_mM,
  "OD"         = od_raw,
  "V (mM/min)" = round(V_mM_min, 6),
  check.names  = FALSE
))

# ----- Kinetic Parameters from Lineweaver-Burk -----
Km   <- 1.3144   # mM  — from LB x-intercept
Vmax <- 0.0015   # mM/min — from LB y-intercept

cat(sprintf("\nUsing LB-derived parameters: Km = %.4f mM | Vmax = %.6f mM/min\n", Km, Vmax))

# ============================================================
# PLOT 1 — MICHAELIS-MENTEN CURVE
# ============================================================

# Smooth theoretical MM curve
S_seq  <- seq(0, max(S_mM) * 1.3, length.out = 500)
V_theo <- (Vmax * S_seq) / (Km + S_seq)
mm_curve_df <- data.frame(S = S_seq, V = V_theo)

# Vmax and Km annotation lines
half_Vmax <- Vmax / 2

mm_plot <- ggplot() +
  
  # Theoretical MM curve
  geom_line(data = mm_curve_df,
            aes(x = S, y = V),
            color = "#E05A2B", linewidth = 1.2) +
  
  # Vmax asymptote (dashed)
  geom_hline(yintercept = Vmax,
             linetype = "dashed", color = "grey50", linewidth = 0.7) +
  
  # Half-Vmax horizontal line
  geom_segment(aes(x = 0, xend = Km,
                   y = half_Vmax, yend = half_Vmax),
               linetype = "dotted", color = "#2B4BA8", linewidth = 0.8) +
  
  # Km vertical drop line
  geom_segment(aes(x = Km, xend = Km,
                   y = 0,   yend = half_Vmax),
               linetype = "dotted", color = "#2B4BA8", linewidth = 0.8) +
  
  # Experimental data points
  geom_point(data = df,
             aes(x = S, y = V),
             color = "#2B4BA8", size = 4, shape = 21,
             fill = "#5B8DEF", stroke = 1.2) +
  
  # Vmax annotation
  annotate("text",
           x = max(S_mM) * 1.25, y = Vmax + 0.00003,
           label = paste0("Vmax = ", round(Vmax, 6), " mM/min"),
           hjust = 1, color = "grey40", size = 3.8, fontface = "italic") +
  
  # Km annotation
  annotate("text",
           x = Km + 0.03, y = half_Vmax * 0.4,
           label = paste0("Km = ", round(Km, 4), " mM"),
           hjust = 0, color = "#2B4BA8", size = 3.8, fontface = "italic") +
  
  # Half-Vmax annotation
  annotate("text",
           x = 0.02, y = half_Vmax + 0.00003,
           label = "Vmax / 2",
           hjust = 0, color = "#2B4BA8", size = 3.5, fontface = "italic") +
  
  scale_x_continuous(breaks = seq(0, 1.4, by = 0.2),
                     limits = c(0, 1.4)) +
  scale_y_continuous(limits = c(0, Vmax * 1.15)) +
  
  labs(
    title    = "Michaelis-Menten Curve — Alkaline Phosphatase",
    subtitle = paste0("Substrate: pNPP  |  ε = 18,800 M⁻¹cm⁻¹  |  t = ",
                      reaction_time, " min  |  Km = ", round(Km, 4),
                      " mM  |  Vmax = ", round(Vmax, 6), " mM/min"),
    x        = "[S]  (mM)",
    y        = "V  (mM / min)"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey50", size = 9),
    axis.title       = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey92")
  )

print(mm_plot)

# ============================================================
# PLOT 2 — EADIE-HOFSTEE PLOT
# V vs V/[S]  |  slope = -Km  |  y-intercept = Vmax
# ============================================================

# Calculate V/[S] for each data point
df$V_over_S <- df$V / df$S

# Linear regression: V = Vmax + (-Km) × (V/[S])
eh_fit        <- lm(V ~ V_over_S, data = df)
eh_slope      <- coef(eh_fit)[2]        # should be ≈ -Km
eh_intercept  <- coef(eh_fit)[1]        # should be ≈ Vmax
eh_r2         <- summary(eh_fit)$r.squared

Km_EH   <- -eh_slope
Vmax_EH <- eh_intercept

cat("\n=== Eadie-Hofstee Parameters ===\n")
cat(sprintf("Slope (−Km)  : %.6f  →  Km   = %.4f mM\n",  eh_slope,     Km_EH))
cat(sprintf("Y-intercept  : %.6f  →  Vmax = %.6f mM/min\n", eh_intercept, Vmax_EH))
cat(sprintf("R²           : %.4f\n", eh_r2))

# Regression line for EH plot
x_eh_seq   <- seq(0, max(df$V_over_S) * 1.15, length.out = 300)
eh_line_df <- data.frame(
  V_over_S = x_eh_seq,
  V        = eh_intercept + eh_slope * x_eh_seq
)

eh_plot <- ggplot() +
  
  # Regression line
  geom_line(data = eh_line_df,
            aes(x = V_over_S, y = V),
            color = "#E05A2B", linewidth = 1.2) +
  
  # Data points
  geom_point(data = df,
             aes(x = V_over_S, y = V),
             color = "#2B4BA8", size = 4, shape = 21,
             fill = "#5B8DEF", stroke = 1.2) +
  
  # Y-intercept point (Vmax)
  geom_point(aes(x = 0, y = Vmax_EH),
             color = "#228B22", size = 4, shape = 17) +
  annotate("text",
           x = max(df$V_over_S) * 0.02,
           y = Vmax_EH + 0.000015,
           label = paste0("Vmax = ", round(Vmax_EH, 6), " mM/min"),
           hjust = 0, color = "#228B22", size = 3.8, fontface = "italic") +
  
  # Slope and Km annotation
  annotate("text",
           x = max(df$V_over_S) * 0.45,
           y = max(df$V) * 0.55,
           label = paste0("Slope = −Km = ", round(eh_slope, 4),
                          "\nKm = ", round(Km_EH, 4), " mM",
                          "\nR² = ", round(eh_r2, 4)),
           hjust = 0, color = "grey25", size = 3.8, fontface = "italic") +
  
  scale_x_continuous(limits = c(0, max(df$V_over_S) * 1.2)) +
  scale_y_continuous(limits = c(0, max(df$V) * 1.2)) +
  
  labs(
    title    = "Eadie–Hofstee Plot — Alkaline Phosphatase",
    subtitle = paste0("Substrate: pNPP  |  ε = 18,800 M⁻¹cm⁻¹  |  t = ",
                      reaction_time, " min  |  Km = ", round(Km_EH, 4),
                      " mM  |  Vmax = ", round(Vmax_EH, 6), " mM/min"),
    x        = "V / [S]  (mM/min per mM  =  min⁻¹)",
    y        = "V  (mM / min)"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey50", size = 9),
    axis.title       = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey92")
  )

print(eh_plot)