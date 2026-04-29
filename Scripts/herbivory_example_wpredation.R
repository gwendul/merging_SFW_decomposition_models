# Herbivory analysis:

library(pacman)
p_load(deSolve, FME, tidyverse, yaml)
verbose = F

# ---- Load models ----
source("R/millennial_model_herbnem.R")
source("R/derive_millennial_parms.R")
source("R/init_millennial_state.R")
source("R/plot_ode_out.R")


parms <- yaml::read_yaml("config/common_herbivore.yml")
parms  <- modifyList(parms, yaml::read_yaml("config/millennial.yml"))
parms  <- modifyList(parms, yaml::read_yaml("config/millennial_herbnem.yml"))

# All roots to mineral soil:
parms$root_to_organic = 0 

parms <- derive_millennial_parms(parms)

parms$wood_mortality = 0

# With herbivores:
millennial_herbnem = ode(
  times = seq(1, 365*2000, by = 1),
  y     = init_millennial_state(HerbNem = T),
  func  = millennial_model_herbnem,
  parms = parms
)

plot_ode_output(millennial_herbnem,
                variable_cols = names(init_millennial_state(HerbNem = T)))

head(millennial_herbnem)
tail(millennial_herbnem)

write.csv(millennial_herbnem, "millennial_herbnem_wpred.csv")


# End Timepoint/Equilibrium
times <- seq(1, 365*2000, by = 1)

# looking at the last few years to check for stability
timeCheck_stable <- tibble(data.frame(millennial_herbnem)) %>%
  filter(
    time %in% c(
      max(times),
      max(times)-365,
      max(times)-365*2,
      max(times)-365*3
    )
  ) 

# ending parameters for model
timeEND <- as.data.frame(millennial_herbnem) %>%
  slice_tail(n = 1) %>%
  select(-time) %>%
  unlist()

write.csv(timeEND, "timeEND_stable_herbivore_wpred.csv")
  

# new simulation
times2 <- seq(1, 365*10, by = 1)

new_herbnem <- ode(
  y = timeEND,
  times = times2,
  func = millennial_model_herbnem,
  parms = parms
)

plot_ode_output(new_herbnem,
                variable_cols = names(init_millennial_state(HerbNem = T)))

write.csv(new_herbnem, "new_herbnem_wpred.csv")


# Simulate the scenarios (with and without nematodes)
timeEND2 <- as.data.frame(new_herbnem) %>%
  slice_tail(n = 1) %>%
  select(-time) %>%
  unlist()

write.csv(timeEND, "timeEND2_scenarios_herbivore_wpred.csv")

fullModel <- timeEND2
noPred    <- replace(timeEND2, "Predator", 0) #setting predators = 0 
noNematode    <- replace(timeEND2, c("Predator", "RootHerb"), 0) #setting herbivore = 0 

# Run 10 year scenarios
out_full <- ode(y = fullModel,  
                times = times2, 
                func = millennial_model_herbnem, 
                parms = parms
                )

write.csv(out_full, "out_full_wpred.csv")

out_noNematode <- ode(y = noNematode,
                       times = times2, 
                       func = millennial_model_herbnem, 
                       parms = parms
                       )

write.csv(out_noNematode, "out_noNematode.csv")

out_noPred <- ode(y = noPred, 
                  times = times2, 
                  func = millennial_model_herbnem, 
                  parms = parms
                  )

write.csv(out_noPred, "out_noPred.csv")


# Run the three scenarios - 100 years
times3 <- seq(1, 365*100, by = 1)

out_full100 <- ode(y = fullModel,  
                   times = times3, 
                   func = millennial_model_herbnem, 
                   parms = parms
                   )

write.csv(out_full100, "out_full_wpred_100.csv")

out_noNematode100 <- ode(y = noNematode,
                          times = times3,
                          func = millennial_model_herbnem, 
                          parms = parms
                          )

write.csv(out_noNematode100, "out_noNematode_100.csv")

out_noPred100 <- ode(y = noPred, 
                  times = times3, 
                  func = millennial_model_herbnem, 
                  parms = parms
)

write.csv(out_noPred100, "out_noPred_100.csv")


vars <- colnames(out_full)[2:15]


# Helper to reshape an ode output the same way the function does internally
to_long <- function(out, scenario) {
  df <- as.data.frame(out)
  names(df)[1] <- "time"
  df |>
    select(time, all_of(vars)) |>
    pivot_longer(-time, names_to = "state", values_to = "value") |>
    mutate(scenario = scenario)
}

# Base plot from function (full food web)
p <- plot_ode_output(out_full, variable_cols = vars)

# other two scenarios
p <- p +
  geom_line(data = to_long(out_full, "Full Food Web"), 
            aes(time, value, color = scenario)) +
  geom_line(data = to_long(out_noPred, "No Predators"),
            aes(time, value, color = scenario)) +
  geom_line(data = to_long(out_noNematode, "No Nematodes"),
            aes(time, value, color = scenario)) +
  scale_color_manual(values = c("Full Food Web" = "black",
                                "No Predators"   = "blue",
                                "No Nematodes" = "lightgreen")) +
  labs(color = "Scenario")
p
ggsave("simulations_fullfoodweb.png",  plot = p,  width = 9.5, height = 6, dpi = 300)

#plot the 100 year simulation scenarios, but 
# Keep only one day per year (default: day 1)
keep_yearly <- function(out, day_of_year = 1) {
  out[out[, "time"] %% 365 == day_of_year, ]
}

p2 <- plot_ode_output(keep_yearly(out_full100), variable_cols = vars)

p2 <- p2 +
  geom_line(data = to_long(keep_yearly(out_full100), "Full Food Web"), 
            aes(time, value, color = scenario)) +
  geom_line(data = to_long(keep_yearly(out_noPred100), "No Predators"),
            aes(time, value, color = scenario)) +
  geom_line(data = to_long(keep_yearly(out_noNematode100), "No Nematodes"),
            aes(time, value, color = scenario)) +
  scale_color_manual(values = c("Full Food Web" = "black",
                                "No Predators"   = "blue",
                                "No Nematodes" = "lightgreen")) +
  labs(color = "Scenario")
p2

ggsave("simulations_fullfoodweb_100yrs.png",  plot = p2,  width = 9.5, height = 6, dpi = 300)


# Combine the three scenario outputs into one long data frame
to_long <- function(out, scenario) {
  as.data.frame(out) |>
    pivot_longer(-time, names_to = "Pool", values_to = "Value") |>
    mutate(Scenario = scenario)
}

scenario_long <- bind_rows(
  to_long(out_full,      "Full Food Web"),
  to_long(out_noNematode,    "No Nematodes"),
  to_long(out_noPred, "No Predators")
)

scenario_colors <- c("Full Food Web" = "black",
                     "No Predators"   = "blue",
                     "No Nematodes" = "lightgreen")


# ---- Fig 1: endpoint (max time) difference from Full model ----
full_endpoint <- scenario_long %>%
  filter(Scenario == "Full Food Web", time == max(time)) %>%
  select(Pool, full_val = Value)

diff_end <- scenario_long %>%
  filter(Scenario != "Full Food Web", time == max(time)) %>%
  left_join(full_endpoint, by = "Pool") %>%
  mutate(Diff = Value - full_val)

fig_end <- ggplot(diff_end, aes(x = Scenario, y = Diff, fill = Scenario)) +
  geom_col(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  facet_wrap(~ Pool, scales = "free_y") +
  scale_fill_manual(values = scenario_colors) +
  labs(y = expression(Delta ~ "endpoint pool size"),
       x = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

fig_end
ggsave("fig_end_fullfoodweb.png",  plot = fig_end,  width = 8, height = 6, dpi = 300)

# ---- Fig 2: whole-trajectory mean difference from Full model ----
full_mean <- scenario_long %>%
  filter(Scenario == "Full Food Web") %>%
  group_by(Pool) %>%
  summarise(full_mean = mean(Value), .groups = "drop")

diff_mean <- scenario_long %>%
  filter(Scenario != "Full Food Web") %>%
  group_by(Scenario, Pool) %>%
  summarise(mean_value = mean(Value), .groups = "drop") %>%
  left_join(full_mean, by = "Pool") %>%
  mutate(Diff = mean_value - full_mean)

fig_mean <- ggplot(diff_mean, aes(x = Scenario, y = Diff, fill = Scenario)) +
  geom_col(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  facet_wrap(~ Pool, scales = "free_y") +
  scale_fill_manual(values = scenario_colors) +
  labs(y = expression(Delta ~ "mean pool size"),
       x = NULL, fill = NULL) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

fig_mean
ggsave("fig_mean_fullfoodweb.png", plot = fig_mean, width = 8, height = 6, dpi = 300)


#effect sizes

# Per-scenario, per-pool summary (endpoint + time-mean)
summary_df <- scenario_long %>%
  group_by(Scenario, Pool) %>%
  summarise(
    endpoint  = Value[time == max(time)],
    time_mean = mean(Value),
    .groups = "drop"
  )

# Baseline values (Full model)
baseline <- summary_df %>%
  filter(Scenario == "Full Food Web") %>%
  select(Pool, baseline_end = endpoint, baseline_mean = time_mean)

# Effect sizes for the removal scenarios
effect_sizes <- summary_df %>%
  filter(Scenario != "Full Food Web") %>%
  left_join(baseline, by = "Pool") %>%
  mutate(
    # Endpoint effects
    delta_end       = endpoint - baseline_end,
    pct_end         = 100 * delta_end / baseline_end,
    LRR_end         = log(endpoint / baseline_end),
    # Time-averaged effects
    delta_mean      = time_mean - baseline_mean,
    pct_mean        = 100 * delta_mean / baseline_mean,
    LRR_mean        = log(time_mean / baseline_mean)
  ) %>%
  select(Scenario, Pool,
         baseline_end, endpoint, delta_end, pct_end, LRR_end,
         baseline_mean, time_mean, delta_mean, pct_mean, LRR_mean)

effect_sizes
write.csv(effect_sizes, "effect_sizes_fullfoodweb.csv", row.names = FALSE)

effect_sizes <- effect_sizes %>%
  filter(!(Scenario == "No Predators" & Pool %in% c("Predator", "CWD")),
         !(Scenario == "No Nematodes" & Pool %in% c("Predator", "RootHerb", "CWD")))

library(scales)

percent_change <- ggplot(effect_sizes, aes(x = Scenario, y = Pool, fill = pct_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%+.1f%%", pct_mean)), size = 3) +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C",
    midpoint = 0,
    limits = c(-50, 50),     # cap the color scale here
    oob = squish,            # anything beyond limits gets the extreme color
    name = "% change\n(time-mean)"
  ) +
  labs(x = NULL, y = NULL,
       title = "Effect of soil animals on soil C pools") +
  theme_minimal() +
  theme(panel.grid = element_blank())

percent_change

ggsave("percent_change_fullfoodweb.png", plot = percent_change, width = 4.5, height = 7, dpi = 300)
