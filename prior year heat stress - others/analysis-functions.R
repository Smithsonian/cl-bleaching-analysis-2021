##################################################################
##                      Analysis functions                      ##
##                   By: Migdonio A. González                   ##
##                   Edited by: Sean Connolly                   ##
##################################################################
select_weighting_matrix <- function(formula1, 
                                    formula2, 
                                    columns, 
                                    bleach_data, 
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
                bleach_reduced <- bleach_data[, c("bin.score", columns)]
                
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

##################################################################
##                          Variograms                          ##
##################################################################
create_variogram_object_RA <- function(year, bleach_data, residuals){
        object <-  variog(coords = bleach_data[, c("X", "Y")],
                          data = residuals,
                          breaks = seq(0, max(dist(bleach_data[, c("X","Y")])), l = 50))
        
        return(object)
}

##################################################################
##                    Predictions - Response                    ##
##################################################################
generate_prediction_dataframe <- function(column.names, bleach_data, nevs){
        
        pred_range <- seq_range(range(bleach_data[, column.names[1]]), n = 500)
        newdata <- data.frame(DHW = rep(pred_range,4),
                              DHW.X = rep(c(0,3,6,9), each = 500))
        
        colnames(newdata) <- column.names
        
        # Dataframe with zeros.
        m <- matrix(0, ncol = nevs, nrow = 4 * 500)
        m <- data.frame(m)
        colnames(m) <- paste0("X", 1:nevs)
        
        # Unified newdata dataframe. 
        newdata_SAC <- cbind(newdata, m)
        
        # Returning dataframe. 
        return(list(no.ev = newdata, ev = newdata_SAC))
}

generate_predictions <- function(column.names, model, newdata, level = 0){
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
        pr_bleach_df <- data.frame(DHW = newdata[,column.names[1]],
                                   DHWM.X = factor(newdata[,column.names[2]]),
                                   fit = pr_bleach_fit$fit,
                                   response = pr_bleach_response$fit,
                                   se = pr_bleach_se,
                                   lowerCI = pr_bleach_ci.l,
                                   upperCI = pr_bleach_ci.h)
        
        colnames(pr_bleach_df)[1:2] <- column.names
        
        # Return.
        return(pr_bleach_df)
        
}

create_plot_object <- function(column.names, plot.title, legend.title, x.label = " ", pred_table){
        p <- ggplot(data = pred_table, aes_string(x = column.names[1], y = "response", group = column.names[2])) + 
                geom_line(aes_string(color = column.names[2]), size = 1) + 
                scale_x_continuous(breaks = seq(0, 14, 2), limits = c(0,14), expand = c(0,0)) +
                scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
                geom_ribbon(aes_string(ymin = "lowerCI", ymax = "upperCI", fill = column.names[2]), alpha = 0.2) +
                labs(x = "DHW (°C-weeks)", y = "Probability of severe bleaching", title = plot.title, col = legend.title, fill = legend.title) +
                theme(legend.position = c(0,0.8), legend.title = element_text(size = 9),
                      legend.text = element_text(size = 8), legend.key.size = unit(0.9, "lines"),
                      legend.background = element_blank(), legend.justification = "left", 
                      axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank()) +
                scale_fill_manual(values=c("deepskyblue1", "darkgoldenrod1", "red", "black"), aesthetics = c("color", "fill"))
        
        return(p)
}