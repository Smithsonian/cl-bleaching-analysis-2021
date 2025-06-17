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
                bleach_reduced <- bleach_data[, c("bin.score",columns,"region")]
                
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
                                       residuals){
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
                                          nevs = 25){
        # Central
        pred_range <- seq_range(range(bleach_data[bleach_data$region == "Central", columns]), n = 500)
        newdata_central <- data.frame(DHW = pred_range,
                                      region = factor("Central"))
        # North
        pred_range <- seq_range(range(bleach_data[bleach_data$region == "North", columns]), n = 500)
        newdata_north <- data.frame(DHW = pred_range,
                                    region = factor("North"))
        # South
        pred_range <- seq_range(range(bleach_data[bleach_data$region == "South", columns]), n = 500)
        newdata_south <- data.frame(DHW = pred_range,
                                    region = factor("South"))
        
        # Dataframe with zeros.
        m <- matrix(0, ncol = nevs, nrow = 3 * 500)
        m <- data.frame(m)
        colnames(m) <- paste0("X", 1:nevs)
        
        # Unified newdata dataframe. 
        newdata_region <- bind_rows(newdata_central,
                                    newdata_north,
                                    newdata_south)
        colnames(newdata_region)[1] <- columns
        
        newdata_region_SAC <- cbind(newdata_region, m)
        colnames(newdata_region_SAC)[1] <- columns
        
        # Returning dataframe. 
        return(list(no.ev = newdata_region, ev = newdata_region_SAC))
}

generate_predictions <- function(model, 
                                 newdata, 
                                 columns, 
                                 level = 0){
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
                                   region = newdata$region,
                                   fit = pr_bleach_fit$fit,
                                   response = pr_bleach_response$fit,
                                   se = pr_bleach_se,
                                   lowerCI = pr_bleach_ci.l,
                                   upperCI = pr_bleach_ci.h)
        
        colnames(pr_bleach_df)[1] <- columns
        
        # Return.
        return(df = pr_bleach_df)
        
}

#x = "DHW (°C-week)"
#x = " "

create_plot_object <- function(year, 
                               columns, 
                               x.label, 
                               plot.title, 
                               pred_table){
        p <- ggplot(data = pred_table, aes_string(x = columns, y = "response", group = "region")) + 
                geom_line(aes(color = region), size = 1) + 
                scale_x_continuous(breaks = seq(0, 14, 2), limits = c(0,14), expand = c(0,0)) +
                scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0,1), expand = c(0,0)) +
                geom_ribbon(aes(ymin = lowerCI, ymax = upperCI, fill = region), alpha = 0.2) +
                labs(x = x.label, y = "Probability of severe bleaching", title = plot.title, col = "Region", fill = "Region") +
                theme(legend.position = c(0, 0.8), legend.title = element_text(size = 9),
                      legend.text = element_text(size = 8), legend.key.size = unit(0.7, "lines"),
                      legend.background = element_blank(), legend.justification = "left", 
                      axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
                      axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      #panel.border = element_blank(),
                      panel.background = element_blank()) +
                scale_fill_manual(values=c("red", "darkgoldenrod1", "deepskyblue1"), aesthetics = c("color", "fill"), labels = c("Northern", "Central", "Southern"))+
                annotate("text", x = 13, y = 0.125, label = year) #+ theme(legend.position="none")
        
        return(p)
}