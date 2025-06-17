##################################################################
##                      Analysis functions                      ##
##                   By: Migdonio A. González                   ##
##                   Edited by: Sean Connolly                   ##
##################################################################
select_weighting_matrix <- function(bleach_data, 
                                    columns,
                                    formula1,
                                    candidates_list)
{
        # Empty dataframe. 
        summary.df <- data.frame(candidate = numeric(), 
                                 nev = numeric(), 
                                 p = numeric(), 
                                 morans = numeric(),
                                 vector.index = numeric())
        
        # Total candidates in this set. 
        print(length(candidates_list))
        
        # For loop.
        for (i in 1:length(candidates_list)){
                # Bleach data reduced.
                bleach_reduced <- bleach_data[, c("bin.score",columns)]
                
                # Using ME to calculate the eigenvectors.
                ME.object <- ME(formula1, 
                                data = bleach_data,
                                family = binomial,
                                listw = candidates_list[[i]])
                
                bleach_reduced <- cbind(bleach_reduced, ME.object$vectors)
                evs <- ncol(ME.object$vectors)
                
                # Extracting p.value.
                p.value <- unname(tail(ME.object$selection[,3], n  = 1))
                
                # Moran's I statistic.  
                morans.I <- unname(tail(ME.object$selection[,2], n  = 1))
                
                # Extracting vector index.
                index <- paste(unname(ME.object$selection[2:nrow(ME.object$selection),1]), collapse = " ")
                
                # Building the return dataframe.
                dfr <- data.frame(candidate = i,
                                  nev = evs,
                                  p = p.value,
                                  moran = morans.I,
                                  vector.index = index)
                summary.df <- rbind(summary.df, dfr)
                
                # Tracking progress.
                print(i)
                print(paste("The number of eigenvectors for this run is:", evs))
        }
        
        
        # Returning dataframe.
        return(summary.df)
}

create_mem_dataframe <- function(mems, 
                                 subset.year, 
                                 nrow = 4128)
{
        mem.df <- data.frame(matrix(0, ncol = ncol(mems), nrow = nrow))
        mem.df[subset.year, ] <- mems
        return(mem.df)
}

##################################################################
##                          Variograms                          ##
##################################################################
create_variogram_object_RA <- function(year, 
                                       bleach_data, 
                                       residuals)
{
        object <-  variog(coords = bleach_data[, c("X", "Y")],
                          data = residuals,
                          breaks = seq(0, max(dist(bleach_data[, c("X","Y")])), l = 50))
        
        return(object)
}

##################################################################
##                    Predictions - Response                    ##
##################################################################
generate_prediction_dataframe <- function(bleach_data, 
                                          columns, 
                                          nevs)
{
        
        # Year
        pred_range <- seq_range(range(bleach_data[,columns]), n = 500)
        newdata_year <- data.frame(DHW = pred_range)
        colnames(newdata_year) <- columns
        
        # Dataframe with zeros. 
        m <- matrix(0, ncol = nevs, nrow = 500)
        m <- data.frame(m)
        colnames(m) <- paste0("X",seq(1,nevs))
        
        # Dataframe with eigenvectors columns (these are filled with zeros).
        newdata_year_SAC <- cbind(newdata_year, m)
        colnames(newdata_year_SAC)[1] <- columns
        
        # Returning dataframe. 
        return(list(no.ev = newdata_year, ev = newdata_year_SAC))
}

generate_prediction_dataframe_full <- function(bleach_data, 
                                               nevs)
{
        # 1998
        pred_range <- seq_range(range(bleach_data$DHW[bleach_data$year == "1998"]), n = 500)
        newdata_year_1998 <- data.frame(DHW = pred_range,
                                        year = factor("1998"))
        # 2002
        pred_range <- seq_range(range(bleach_data$DHW[bleach_data$year == "2002"]), n = 500)
        newdata_year_2002 <- data.frame(DHW = pred_range,
                                        year = factor("2002"))
        # 2016
        pred_range <- seq_range(range(bleach_data$DHW[bleach_data$year == "2016"]), n = 500)
        newdata_year_2016 <- data.frame(DHW = pred_range,
                                        year = factor("2016"))
        
        # 2017
        pred_range <- seq_range(range(bleach_data$DHW[bleach_data$year == "2017"]), n = 500)
        newdata_year_2017 <- data.frame(DHW = pred_range,
                                        year = factor("2017"))
        
        # 2020
        pred_range <- seq_range(range(bleach_data$DHW[bleach_data$year == "2020"]), n = 500)
        newdata_year_2020 <- data.frame(DHW = pred_range,
                                        year = factor("2020"))
        
        # Dataframe with zeros. 25 zeros because there are 25 eigenvectors.
        m <- matrix(0, ncol = nevs, nrow = 5 * 500)
        m <- data.frame(m)
        colnames(m) <- paste0("X",seq(1,nevs))
        
        # Unified newdata dataframe. 
        newdata_year <- bind_rows(newdata_year_1998,
                                  newdata_year_2002,
                                  newdata_year_2016,
                                  newdata_year_2017,
                                  newdata_year_2020)
        newdata_year_SAC <- cbind(newdata_year, m)
        
        # Returning dataframe. 
        return(list(no.ev = newdata_year, ev = newdata_year_SAC))
}

generate_predictions <- function(model, 
                                 newdata,
                                 columns,
                                 level = 0)
{
        # Fit. Linear Scale.
        pr_bleach_fit <- predict(model, 
                                 newdata = newdata,
                                 level = level,
                                 se.fit = TRUE)
        # Response.
        pr_bleach_response <- predict(model, 
                                      newdata = newdata,
                                      level = level,
                                      type = "response",
                                      se.fit = TRUE)
        # Standard errors. 
        pr_bleach_se <- pr_bleach_fit$se.fit
        
        # Confidence Intervals.
        pr_bleach_ci.l <- inv.logit(pr_bleach_fit$fit - 1.96 * pr_bleach_se)
        
        pr_bleach_ci.h <- inv.logit(pr_bleach_fit$fit + 1.96 * pr_bleach_se)
        
        # Return Dataframe.
        pr_bleach_df <- data.frame(DHW = newdata[,columns],
                                   fit = pr_bleach_fit$fit,
                                   response = pr_bleach_response$fit,
                                   se = pr_bleach_se,
                                   lowerCI = pr_bleach_ci.l,
                                   upperCI = pr_bleach_ci.h)
        colnames(pr_bleach_df)[1] <- columns
        
        # Return.
        return(pr_bleach_df)
        
}

generate_predictions_full <- function(model, 
                                      newdata, 
                                      level = 0)
{
        # Fit. Linear Scale.
        pr_bleach_fit <- predict(model, 
                                 newdata = newdata,
                                 level = level,
                                 se.fit = TRUE)
        # Response.
        pr_bleach_response <- predict(model, 
                                      newdata = newdata,
                                      level = level,
                                      type = "response",
                                      se.fit = TRUE)
        # Standard errors. 
        pr_bleach_se <- pr_bleach_fit$se.fit
        
        # Confidence Intervals.
        pr_bleach_ci.l <- inv.logit(pr_bleach_fit$fit - 1.96 * pr_bleach_se)
        
        pr_bleach_ci.h <- inv.logit(pr_bleach_fit$fit + 1.96 * pr_bleach_se)
        
        # Return Dataframe.
        pr_bleach_df <- data.frame(DHW = newdata$DHW,
                                   year = newdata$year,
                                   fit = pr_bleach_fit$fit,
                                   response = pr_bleach_response$fit,
                                   se = pr_bleach_se,
                                   lowerCI = pr_bleach_ci.l,
                                   upperCI = pr_bleach_ci.h)
        
        # Return.
        return(pr_bleach_df)
        
}

create_plot_object <- function(plot.title, 
                               columns, 
                               x.label = "", 
                               pred_table)
{
        p <- ggplot(data = pred_table, aes_string(x = columns, y = "response")) + 
                geom_line() + 
                scale_x_continuous(breaks = seq(0, 14, 2), limits = c(0,14)) +
                geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha = 0.2) +
                labs(x = x.label, y = "Pr. bleach", title = plot.title) +
                theme(legend.position = c(0.7, 0.75), legend.title = element_text(size = 9),
                      legend.text = element_text(size = 8), legend.key.size = unit(0.9, "lines"),
                      legend.background = element_blank(), legend.justification = "left", 
                      panel.background = element_rect(fill = "white", colour = "black"),
                      panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
                      panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
                scale_fill_manual(values=c("red", "green", "blue","black","orange"), aesthetics = c("color", "fill"))
        
        return(p)
}

#x = "DHW (°C-week)"
#x = " "
create_plot_object_full <- function(plot.title, 
                                    x.label = "", 
                                    pred_table)
{
        p <- ggplot(data = pred_table, aes(x = DHW, y = response, group = year)) + 
                geom_line(aes(color = year), size = 1) + 
                scale_x_continuous(breaks = seq(0, 14, 2), limits = c(0,14), expand = c(0,0)) +
                scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
                geom_ribbon(aes(ymin = lowerCI, ymax = upperCI, fill = year), alpha = 0.2) +
                labs(x = x.label, y = "Probability of severe bleaching", title = plot.title, col = "Year",fill = "Year") +
                theme(legend.position = c(0, 0.8), legend.title = element_text(size = 9),
                      legend.text = element_text(size = 8), legend.key.size = unit(0.9, "lines"),
                      legend.background = element_blank(), legend.justification = "left", 
                      axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      #panel.border = element_blank(),
                      panel.background = element_blank()) +
                scale_fill_manual(values=c("red", "darkgoldenrod1", "black", "deepskyblue1","purple"), aesthetics = c("color", "fill"))
        
        return(p)
}