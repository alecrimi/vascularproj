#####################################################################################################################################################
#####################################################################################################################################################
###                                                                                                                                               ###
###    Created on:        20th March 2019                                                                                                         ###
###    Last modified on:  12th February 2020                                                                                                      ###
###    Lukas Frick                                                                                                                                ###
###                                                                                                                                               ###
###                                                                                                                                               ###
###    **********************************************************************************************************************************         ###
###    ***                                                                                                                            ***         ###
###    ***   Instructions for running the script                                                                                      ***         ###
###    ***                                                                                                                            ***         ###
###    ***   Run all code by selecting 'Code' > 'Source with Echo', or hitting 'Ctrl + Shift + Enter'.                                ***         ###
###    ***   A Windows dialog will appear for selecting a folder, but it will appear in the background!                               ***         ###
###    ***   Select it by clicking on its tab on the Windows taskbar.                                                                 ***         ###
###    ***   Use the dialog to select either a folder containing .raw files, or a folder (or hard drive) containing                   ***         ###
###    ***   one or more folders that contain .raw files.                                                                             ***         ###
###    ***                                                                                                                            ***         ###
###    ***   That's it, the script should run automatically.                                                                          ***         ###
###    ***   You can follow its progress while it's running by reading the messages in the RStudio console,                           ***         ###
###    ***   or by opening the log text file, which is located in the "Stitch_log" folder created in the folder you selected.         ***         ###
###    ***                                                                                                                            ***         ###
###    ***   Please report unexpected errors by sending me an e-mail (lukas.frick@usz.ch).                                            ***         ###
###    ***                                                                                                                            ***         ###
###    **********************************************************************************************************************************         ###
###                                                                                                                                               ###
###                                                                                                                                               ###
###    This R script allows fully automatic stitching of all .raw files in a directory.                                                           ###
###    It requires the knowledge of the layout of tiles (rows and columns),                                                                       ###
###    and also of the order in which the tiles of this layout are acquired by the SPIM                                                           ###
###    (see the section 'Physical z stack tile layout' below).                                                                                    ###
###    It assumes that the last modified timestamp of the .raw files corresponds to the time of acquisition.                                      ###
###    It can extract the information on voxel dimensions and the number of z planes from the metadata .txt files.                                ###
###    It may be possible to extract the mechanical displacement between tiles from the metadata,                                                 ###
###    but in the current version of the script, it is predetermined (i.e., 5 millimeters).                                                       ###
###                                                                                                                                               ###
###    It uses the system shell (the Windows command line) to call external programs, namely Fiji and Terastitcher/Teraconverter.                 ###
###    Standalone (portable) versions of these programs are called.                                                                               ###
###                                                                                                                                               ###
###                                                                                                                                               ###
#####################################################################################################################################################
#####################################################################################################################################################



#####################################################################################################################################################
#####################################################################################################################################################
###                                                                                                                                               ###
###    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                                                     ###
###    !!!                                !!!                                                                                                     ###
###    !!!  Physical z stack tile layout  !!!                                                                                                     ###
###    !!!                                !!!                                                                                                     ###
###    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                                                     ###
###                                                                                                                                               ###
###                                                                                                                                               ###
###         #############                                                                                                                         ###
###         #     #     #                                                                                                                         ###
###         #  5  #  4  #                                                                                                                         ###
###         #     #     #                                                                                                                         ###
###         #############                                                                                                                         ###
###         #     #     #                                                                                                                         ###
###         #  6  #  3  #                                                                                                                         ###
###         #     #     #           "1" is the oldest file, and "8" is the newest file                                                            ###
###         #############            (according to the 'last modified' timestamp)                                                                 ###
###         #     #     #                                                                                                                         ###
###         #  7  #  2  #                                                                                                                         ###
###         #     #     #                                                                                                                         ###
###         #############                                                                                                                         ###
###         #     #     #                                                                                                                         ###
###         #  8  #  1  #                                                                                                                         ###
###         #     #     #                                                                                                                         ###
###         #############                                                                                                                         ###
###                                                                                                                                               ###
#####################################################################################################################################################
#####################################################################################################################################################






# Load packages -----------------------------------------------------------

library("processx")




# Define paths ------------------------------------------------------------

ImageJ_app_directory  <- "~/Automatic_stitch/fiji-win64"







# Define global lookup vectors --------------------------------------------

physical_layout_mat <- cbind(5:8, 4:1)


emission_filter_search_strings <- c(
  "498"   = "_498_509-22_",
  "515LP" = "_515LP_",
  "565"   = "_565_585-40_",
  "405"   = "_405-488-647-Tripleblock_2x_",
  "405D"  = "_405_",
  "594"   = "_594_",
  "647"   = "_647_nm_405-488-561-640-Quadrupleblock_",
  "520SG" = "_520SG_",
  "562SG" = "_562SG_"
)

sample_rotation_search_strings <- c(
  "0"   = "_(rot_0|0deg|0_deg)_",
  "180" = "_(rot_180|180deg|180_deg)_"
)

known_values <- list(
  "pixel_size_in_microns"  = c(3.26),
  "z_step_size_in_microns" = c(3),
  "number_of_z_planes"     = 2500L,
  "laser_power_in_percent" = c(8, 10),
  "laser_wavelength"       = "488 nm"
)

known_file_sizes <- c(3.906, 0.516, 0.523, 1.297, 1.82, 2.078, 2.602, 3.383, 3.641, 3.516)


### logging ###
log_file_path <- NULL
log_time_zero <- NULL

### Terastitcher logging ###
captured_output_list <- list()
printed_output_list <- list()
current_output <- c()



# General helper functions ------------------------------------------------

TidyPaths <- function(paths_vec) {
  results_vec <- gsub("\\/", "/", paths_vec,   fixed = TRUE)
  results_vec <- gsub("//", "/",  results_vec, fixed = TRUE)
  results_vec <- gsub("\\", "/",  results_vec, fixed = TRUE)
  return(results_vec)
}


dhm <- function(time_in_seconds) {
  t <- round(time_in_seconds / 60)
  num_days <- t %/% (60L * 24L)
  num_hours <- t %/% 60L %% 24L
  num_minutes <- t %% 60L
  result_vec <- paste0(paste0(formatC(num_hours, width = 2, format = "d", flag = "0"), "h"),
                       " ",
                       paste0(formatC(num_minutes, width = 2, format = "d", flag = "0"), "min")
                       )
  result_vec <- ifelse(num_days > 0, paste0(num_days, "d ", result_vec), result_vec)
  return(result_vec)
}

dhms <- function(t, drop_hours = TRUE, add_space = FALSE) {
  num_days <- t %/% (60L * 60L * 24L)
  num_hours <- t %/% (60L * 60L) %% 24
  num_minutes <- t %/% 60L %% 60
  num_seconds <- t %% 60L
  result_vec <- paste0(paste0(formatC(num_minutes, width = 2, format = "d", flag = "0"), if (add_space) " min" else "min"),
                       " ",
                       paste0(formatC(num_seconds, width = 2, format = "d", flag = "0"), if (add_space) " sec" else "s")
                       )
  hours_added_vec <- paste0(num_hours, if (add_space) " h " else "h ", result_vec)
  if (drop_hours) {
    result_vec <- ifelse((num_hours > 0) | (num_days > 0), hours_added_vec, result_vec)
  } else {
    result_vec <- hours_added_vec
  }
  result_vec <- ifelse(num_days > 0, paste0(num_days, if (add_space) " d " else "d ", result_vec), result_vec)
  return(result_vec)
}







# Functions for locating data folders -------------------------------------

ContainsRawFiles <- function(check_path, raw_files_regex = "\\.raw$") {
  stopifnot(dir.exists(check_path))
  any(grepl(raw_files_regex, list.files(check_path)))
}


FindFoldersWithRawFiles <- function(root_directory) {
  all_folder_paths <- list.dirs(root_directory, recursive = FALSE, full.names = TRUE)
  contain_raw_files <- vapply(all_folder_paths, ContainsRawFiles, logical(1))
  results_vec <- all_folder_paths[contain_raw_files]
  results_vec <- TidyPaths(results_vec)
  return(results_vec)
}







# Functions for creating the hierarchy of folders -------------------------

FormatFolderNumber <- function(number_of_folders, displacement_in_microns) {
  formatC((seq_len(number_of_folders) - 1L) * displacement_in_microns * 10L, flag = "0", width = 6L, format = "d")
}

CreateTerastitcherFolders <- function(num_rows, num_columns, displacement_in_microns, use_directory, create_folders = TRUE) {
  column_names <- FormatFolderNumber(num_columns, displacement_in_microns)
  row_names <- FormatFolderNumber(num_rows, displacement_in_microns)
  lowest_folder_paths <- c()
  for (column_name in column_names) {
    dir.create(file.path(use_directory, column_name), showWarnings = FALSE)
    row_paths <- file.path(use_directory, column_name, paste0(column_name, "_", row_names))
    lowest_folder_paths <- c(lowest_folder_paths, row_paths)
    if (create_folders) {
      for (row_path in row_paths) {
        dir.create(row_path, showWarnings = FALSE)
      }
    }
  }
  return(lowest_folder_paths)
}


MakeCombinationsMat <- function(emission_filters, rotations) {
  hemispheres_mat <- as.matrix(expand.grid(as.character(rotations), emission_filters, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)[, 2:1])
  colnames(hemispheres_mat) <- c("Emission_filter", "Sample_rotation")
  return(hemispheres_mat)
}


PrepareForSplitting <- function(raw_files_folder_path,
                                split_and_stitch_root_path = NULL,
                                displacement_in_microns    = 5000,
                                use_physical_layout_mat    = physical_layout_mat,
                                emission_filters           = NULL,
                                stitch_folder_prefix       = "R_STITCH_",
                                raw_files_regex            = "\\.raw$",
                                show_excluded_message      = TRUE,
                                create_folders             = TRUE,
                                no_rotation                = FALSE
                                ) {

  if (is.null(split_and_stitch_root_path)) {
    split_and_stitch_root_path <- dirname(raw_files_folder_path)
  }

  if (is.null(emission_filters)) {
    emission_filters <- c("498", "565")
  }

  found_emission_filters <- emission_filters %in% names(emission_filter_search_strings)
  if (!(all(found_emission_filters))) {
    LogAndStop(paste("The following emission filter name", if (sum(!(found_emission_filters)) == 1) " was" else "s were", " not recognized: ",
                     paste0(emission_filters[!(found_emission_filters)], collapse = ", ")
                     )
               )
  }
  filter_search_strings <- emission_filter_search_strings[emission_filters]

  all_files <- list.files(raw_files_folder_path)
  all_files <- all_files[order(file.mtime(file.path(raw_files_folder_path, all_files)))]

  are_raw_files <- grepl(raw_files_regex, all_files)

  if (!(any(are_raw_files))) {
    LogAndStop("No .raw files were found in the directory specified!")
  }

  if (no_rotation) {
    rotations <- 0L
  } else {
    rotations <- c(0L, 180L)
  }

  num_rows <- nrow(use_physical_layout_mat)
  num_columns <- ncol(use_physical_layout_mat)
  num_files <- num_rows * num_columns
  file_order <- as.vector(use_physical_layout_mat)

  new_top_folder_name <- paste0(stitch_folder_prefix, gsub(" ", "_", basename(raw_files_folder_path), fixed = TRUE))
  new_top_path <- TidyPaths(file.path(split_and_stitch_root_path, new_top_folder_name))

  hemispheres_mat <- MakeCombinationsMat(emission_filters, rotations)

  assign("delete_hemispheres_mat", hemispheres_mat, envir = globalenv())
  assign("delete_filter_search_strings", filter_search_strings, envir = globalenv())
  assign("delete_all_files", all_files, envir = globalenv())
  assign("delete_sample_rotation_search_strings", sample_rotation_search_strings, envir = globalenv())
  assign("delete_no_rotation", no_rotation, envir = globalenv())
  assign("delete_are_raw_files", are_raw_files, envir = globalenv())

  hemisphere_mat_list <- lapply(seq_len(nrow(hemispheres_mat)), function(x) {

    filter <- hemispheres_mat[x, "Emission_filter"]
    rotation <- hemispheres_mat[x, "Sample_rotation"]

    are_this_filter <- grepl(filter_search_strings[[filter]], all_files, fixed = TRUE)
    if (no_rotation) {
      are_this_rotation <- rep(TRUE, length(all_files))
    } else {
      are_this_rotation <- grepl(sample_rotation_search_strings[[as.character(rotation)]], all_files)
    }

    these_file_names <- all_files[are_raw_files & are_this_filter & are_this_rotation]

    top_path_created <- FALSE

    if (length(these_file_names) != num_files) {
      if (show_excluded_message) {
        MessageAndLog(paste0("For the ", filter, " nm emission filter, and with the sample in the ",
                             rotation, "° position, the expected ", num_files, " .raw files were not found! This condition was left out."
                             )
                      )
      }
      return(NULL)
    } else {
      filter_folder_path <- file.path(new_top_path, paste0("filter_", filter))
      rotation_folder_path <- file.path(filter_folder_path, paste0("rot_", rotation))

      if (create_folders) {
        dir.create(new_top_path,         showWarnings = FALSE)
        dir.create(filter_folder_path,   showWarnings = FALSE)
        dir.create(rotation_folder_path, showWarnings = FALSE)
      }

      lowest_folder_paths <- CreateTerastitcherFolders(num_rows, num_columns, displacement_in_microns, rotation_folder_path, create_folders = create_folders)

      results_mat <- cbind("Old_folder_path" = TidyPaths(raw_files_folder_path),
                           "New_top_path"    = new_top_path,
                           "Emission_filter" = filter,
                           "Sample_rotation" = rotation,
                           "Old_file_name"   = these_file_names[file_order],
                           "New_folder_path" = lowest_folder_paths
                           )
      return(results_mat)

    }
  })

  hemisphere_mat_list <- Filter(Negate(is.null), hemisphere_mat_list)

  if (length(hemisphere_mat_list) == 0) {
    LogAndStop(paste0("In '", raw_files_folder_path, "', no collection of ", num_files, " .raw files was found that met the filtering criteria!"))
  }

  results_mat <- do.call(rbind, hemisphere_mat_list)

  return(results_mat)
}





# Functions for checking input folders for completeness -------------------

CheckAllForCompleteness <- function(raw_files_root_or_folder,
                                    split_and_stitch_root_path = NULL,
                                    emission_filters           = NULL,
                                    stitch_folder_prefix       = "R_STITCH_",
                                    raw_files_regex            = "\\.raw$",
                                    raise_error                = TRUE,
                                    raise_metadata_error       = raise_error,
                                    use_metadata               = TRUE,
                                    no_rotation                = FALSE,
                                    use_physical_layout_mat    = physical_layout_mat
                                    ) {
  if (is.null(emission_filters)) {
    emission_filters <- c("498", "565")
  }
  rotations <- c("0", "180")

  if (ContainsRawFiles(raw_files_root_or_folder)) {
    folder_paths_vec <- raw_files_root_or_folder
    are_subfolders <- FALSE
  } else {
    folder_paths_vec <- FindFoldersWithRawFiles(raw_files_root_or_folder)
    if (length(folder_paths_vec) == 0) {
      LogAndStop("No .raw files were found in the chosen folder, nor within any of its subfolders!")
    }
    are_subfolders <- TRUE
  }
  folder_paths_vec <- TidyPaths(folder_paths_vec)

  hemisphere_df_list <- lapply(folder_paths_vec, function(folder_path) {

    arguments_list <- list(
      folder_path,
      NULL,
      emission_filters        = emission_filters,
      stitch_folder_prefix    = stitch_folder_prefix,
      raw_files_regex         = raw_files_regex,
      create_folders          = FALSE,
      show_excluded_message   = FALSE,
      no_rotation             = no_rotation,
      use_physical_layout_mat = use_physical_layout_mat
    )
    assign("delete_foo_arguments_list", arguments_list, envir = globalenv())

    stitch_mat <- PrepareForSplitting(folder_path,
                                      NULL,
                                      emission_filters        = emission_filters,
                                      stitch_folder_prefix    = stitch_folder_prefix,
                                      raw_files_regex         = raw_files_regex,
                                      create_folders          = FALSE,
                                      show_excluded_message   = FALSE,
                                      no_rotation             = no_rotation,
                                      use_physical_layout_mat = use_physical_layout_mat
                                      )

    assign("delete_stitch_mat", stitch_mat, envir = globalenv())

    stitch_df <- MetaDataForStitchMat(stitch_mat, raise_error = raise_metadata_error, tsv_path = MetaDataTSVPath(folder_path, split_and_stitch_root_path))

    stitch_mat <- unique(stitch_mat[, c("Old_folder_path", "New_top_path", "Emission_filter", "Sample_rotation")])

    complete_mat <- MakeCombinationsMat(emission_filters, rotations)
    stitch_vec <- paste0(stitch_mat[, "Emission_filter"], "__", stitch_mat[, "Sample_rotation"])
    complete_vec <- paste0(complete_mat[, "Emission_filter"], "__", complete_mat[, "Sample_rotation"])
    results_df <- data.frame("Folder_path" = folder_path,
                             complete_mat,
                             "Was_present" = complete_vec %in% stitch_vec,
                             stringsAsFactors = FALSE
                             )
    return(results_df)
  })
  hemisphere_df <- do.call(rbind.data.frame, c(hemisphere_df_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  num_filters <- length(emission_filters)
  MessageAndLog(paste0("Checking data for completeness (i.e., looking for 8 .raw files for ",
                       if (num_filters > 1) paste0("each of the ", num_filters, " emission filters") else "a single emission filter",
                       " and ",
                       length(rotations), " sample rotations, i.e. ",
                       length(emission_filters) * length(rotations) * 8,
                       " in total)..."
                       ))

  if (nrow(hemisphere_df) == 0) {
    LogAndStop(paste("No complete set of 8 files, corresponding to one hemisphere, was found in ",
                     if (are_subfolders) "the subfolders of " else "",
                     "'", raw_files_root_or_folder, "'!"
                     ))
  } else {
    are_present <- hemisphere_df[, "Was_present"]
    are_not_present <- !(are_present)
    if (any(are_present)) {
      unique_folders <- unique(hemisphere_df[are_present, "Folder_path"])
      for (unique_folder in unique_folders) {
        are_this_folder <- hemisphere_df[, "Folder_path"] == unique_folder
        num_present <- sum(are_this_folder & are_present)
        MessageAndLog(paste0("In '" , unique_folder, "', complete data were found for ",
                             num_present, " hemisphere",
                             if (num_present == 1) "" else "s", "."
                             )
                      )
      }
    }
    if (any(are_not_present)) {
      unique_folders <- unique(hemisphere_df[are_not_present, "Folder_path"])
      for (unique_folder in unique_folders) {
        are_this_folder <- hemisphere_df[, "Folder_path"] == unique_folder
        indices <-  which(are_this_folder & are_not_present)
        MessageAndLog(paste0("In '" , unique_folder, "', data was missing for the following condition",
                             if (length(indices) == 1) "" else "s", ":"
                             )
                      )
        for (index in indices) {
          MessageAndLog(paste0("- The ", hemisphere_df[index, "Emission_filter"], " nm emission filter and a sample rotation of ",
                               hemisphere_df[index, "Sample_rotation"], "°."
                               )
                        )
        }
        LogBlankLines(1)
      }
      if (raise_error) {
        LogAndStop("An error was raised due to the incomplete data (as described above).")
      }
    }
  }
  LogBlankLines(2)
  return(invisible(NULL))

}





# Functions for reading metadata ------------------------------------------

ReadMetaData <- function(folder_path, full_file_name) {
  if (!(grepl("\\.raw$", full_file_name))) {
    LogAndStop(paste0("The file name '", full_file_name, "' provided to the ReadMetaData function was invalid - it did not end with '.raw'!"))
  }
  meta_data_path <- TidyPaths(file.path(folder_path, paste0(full_file_name, "_meta.txt")))
  meta_data_vec <- read.table(meta_data_path, header = FALSE, quote = "", sep = "\t", stringsAsFactors = FALSE, blank.lines.skip = FALSE)[, 1]
  are_valid_meta_data <- (meta_data_vec == "") | grepl("^\\[.+] ", meta_data_vec)
  if (!(all(are_valid_meta_data))) {
    LogAndStop(paste0("The metadata passed to the RetrieveMetaItem function contained unexpected values: ",
                      paste0(meta_data_vec[!(are_valid_meta_data)], collapse = ", ")
                      ))
  }
  results_list <- list(
    "file_path"              = RetrieveMetaItem("Metadata for file", meta_data_vec),
    "pixel_size_in_microns"  = as.numeric(RetrieveMetaItem("Pixelsize in um", meta_data_vec)),
    "z_step_size_in_microns" = as.numeric(RetrieveMetaItem("z_stepsize",      meta_data_vec)),
    "number_of_z_planes"     = as.integer(RetrieveMetaItem("z_planes",        meta_data_vec)),
    "laser_wavelength"       = RetrieveMetaItem("Laser", meta_data_vec),
    "laser_power_in_percent" = as.numeric(RetrieveMetaItem("Intensity (%)",   meta_data_vec)),
    "filter"                 = RetrieveMetaItem("Filter", meta_data_vec),
    "zoom"                   = RetrieveMetaItem("Zoom", meta_data_vec),
    "shutter"                = RetrieveMetaItem("Shutter", meta_data_vec)
  )
  return(results_list)
}


EscapeForRegex <- function(string_vec) {
  ### Copied from Hmisc:::escapeRegex
  gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", string_vec)
}


RetrieveMetaItem <- function(item, meta_data_string_vec) {
  this_item_regex <- paste0("^\\[", EscapeForRegex(item), "] ")
  result_indices <- grep(this_item_regex, meta_data_string_vec)
  if (length(result_indices) == 0) {
    assign("debug_meta_data_string_vec", meta_data_string_vec, envir = globalenv())
    LogAndStop("The item '", item, "' was not found, so the RetrieveMetaItem function raised an error.")
  }
  my_result <- meta_data_string_vec[[result_indices[[length(result_indices)]]]]
  my_result <- sub(this_item_regex, "", my_result)
  return(my_result)
}


GetMetaDataForStitchMat <- function(stitch_mat) {
  meta_data_list <- lapply(seq_len(nrow(stitch_mat)),
                           function(x) ReadMetaData(folder_path = stitch_mat[x, "Old_folder_path"], full_file_name = stitch_mat[x, "Old_file_name"])
                           )
  meta_data_df <- do.call(rbind.data.frame, c(meta_data_list, list(stringsAsFactors = FALSE, make.row.names = FALSE)))
  results_df <- data.frame(stitch_mat, meta_data_df, stringsAsFactors = FALSE)
  return(results_df)
}



MetaDataForStitchMat <- function(stitch_mat,
                                 check_metadata         = TRUE,
                                 check_file_names       = TRUE,
                                 check_for_known_values = TRUE,
                                 raise_error            = TRUE,
                                 tsv_path               = NULL,
                                 always_make_tsv        = TRUE
                                 ) {
  meta_df <- GetMetaDataForStitchMat(stitch_mat)
  if (check_metadata) {
    error_found <- FALSE
    AnnounceMetaCheck <- function() {
      MessageAndLog(paste0("The metadata of the files in the '", stitch_mat[1, "Old_folder_path"], "' folder will be checked..."))
    }
    for (variable in names(known_values)) {
      unique_values <- unique(meta_df[, variable])
      if (length(unique_values) != 1) {
        error_found <- TRUE
        AnnounceMetaCheck()
        MessageAndLog(paste0("Inconsistent values for '", variable, "' were found in the .raw files within '",
                             stitch_mat[1, "Old_folder_path"], "': ", paste0(unique_values, collapse = ", ")
                             )
                      )
      }
    }
    if (check_file_names) {
      file_names_from_meta_data <- vapply(strsplit(meta_df[, "file_path"], "/", fixed = TRUE), function(x) x[[length(x)]], "")
      are_identical <- vapply(seq_len(nrow(meta_df)), function(x) identical(file_names_from_meta_data[[x]], as.character(meta_df[x, "Old_file_name"])), logical(1))
      if (!(all(are_identical))) {
        if (!(error_found)) {
          error_found <- TRUE
          AnnounceMetaCheck()
        }
        MessageAndLog("In the following cases, the name of the .raw file and the name recorded in the metadata file were inconsistent:")
        max_num_to_report <- 3
        for (i in which(!(are_identical))[seq_len(min(c(sum(!(are_identical)), max_num_to_report)))]) {
          MessageAndLog(paste0(".raw file name: '",  meta_df[i, "Old_file_name"], "', and name recorded in the metadata: '", file_names_from_meta_data[[i]], "'"))
        }
        if (sum(!(are_identical)) > max_num_to_report) {
          MessageAndLog(paste0("... additional inconsistencies were found, but were omitted for brevity. In total, ", sum(!(are_identical)), " file names were not identical."))
        }
      }
    }
    if (check_for_known_values) {
      for (variable in names(known_values)) {
        values_vec <- meta_df[, variable]
        unexpected_values <- unique(values_vec[!(values_vec %in% known_values[[variable]])])
        if (length(unexpected_values) > 0) {
          if (!(error_found)) {
            error_found <- TRUE
            AnnounceMetaCheck()
          }
          MessageAndLog(paste0("Unexpected values for '", variable, "' were found in the .raw files within '",
                               stitch_mat[1, "Old_folder_path"], "': ", paste0(unexpected_values, collapse = ", ")
                               )
                        )
        }
      }
    }
    if (!(is.null(tsv_path)) && (always_make_tsv || error_found)) {
      formatted_meta_df <- FormatMetaDf(meta_df)
      write.table(formatted_meta_df, file = tsv_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
    if (raise_error && error_found) {
      assign("debug_meta_df", meta_df, envir = globalenv())
      LogAndStop("Due to irregularities in the metadata (as detailed above), an error was raised by the MetaDataForStitchMat function.")
    }
    if (error_found) {
      LogBlankLines(1)
    } else {
      MessageAndLog(paste0("The metadata of the files in the '", stitch_mat[1, "Old_folder_path"], "' folder were checked, and no irregularities were found."))
    }
  }
  return(meta_df)
}



FormatMetaDf <- function(meta_df) {
  metadata_filename_splits <- strsplit(meta_df[, "file_path"], "/", fixed = TRUE)
  metadata_root_paths_vec <- sapply(metadata_filename_splits, function(x) {
    path_length <- length(x)
    if (path_length <= 2) {
      return("")
    } else {
      return(paste0(x[seq_len(path_length - 2)], collapse = "/"))
    }
  })
  file_info_df <- file.info(file.path(meta_df[, "Old_folder_path"], meta_df[, "Old_file_name"]), extra_cols = FALSE)
  results_df <- data.frame(
    meta_df[, c("Emission_filter", "Sample_rotation")],
    "Metadata_root_path"   = metadata_root_paths_vec,
    "Metadata_folder_name" = vapply(metadata_filename_splits, function(x) x[[length(x) - 1]], ""),
    "Actual_folder_name"   = vapply(strsplit(meta_df[, "Old_folder_path"], "/", fixed = TRUE),
                                    function(x) x[[length(x)]],
                                    ""
                                    ),
    "Metadata_file_name"   = vapply(metadata_filename_splits, function(x) x[[length(x)]], ""),
    "Actual_file_name"     = meta_df[, "Old_file_name"],
    "Last_modified"        = as.character(file.mtime(file.path(meta_df[, "Old_folder_path"], meta_df[, "Old_file_name"]))),
    "Scan_duration"        = dhms(as.double(file_info_df[, "mtime"] - file_info_df[, "ctime"], units = "secs")),
    "Time_since_start"     = dhms(as.double(file_info_df[, "mtime"] - min(file_info_df[, "ctime"]), units = "secs"), drop_hours = FALSE),
    "Pixel_size"           = meta_df[, "pixel_size_in_microns"],
    "Z_step_size"          = meta_df[, "z_step_size_in_microns"],
    "Number_of_z_planes"   = meta_df[, "number_of_z_planes"],
    "Laser"                = meta_df[, "laser_wavelength"],
    "Intensity"            = paste0(meta_df[, "laser_power_in_percent"], "%"),
    "Zoom"                 = meta_df[, "zoom"],
    "Shutter"              = meta_df[, "shutter"],
    "Metadata_filter"      = meta_df[, "filter"],
    stringsAsFactors = FALSE
  )
}


MetaDataTSVPath <- function(raw_files_folder_path, split_and_stitch_root_path = NULL, sub_folder = "Stitch_log") {
  if (is.null(split_and_stitch_root_path)) {
    split_and_stitch_root_path <- dirname(raw_files_folder_path)
  }
  file_name <- paste0("metadata_", basename(raw_files_folder_path), ".tsv")
  if (is.null(sub_folder)) {
    path_vec <- c(split_and_stitch_root_path, file_name)
  } else {
    dir.create(TidyPaths(file.path(split_and_stitch_root_path, sub_folder)), showWarnings = FALSE)
    path_vec <- c(split_and_stitch_root_path, sub_folder, file_name)
  }
  result_string <- do.call(file.path, as.list(path_vec))
  result_string <- TidyPaths(result_string)
  return(result_string)
}





# Helper functions for invoking system shell calls ------------------------

FormatPathForShell <- function(path) {
  paste0('"', gsub("/", "\\", path, fixed = TRUE), '"')
}






# Functions for splitting .raw files into .tiffs (with ImageJ) ------------

EscapePathForImageJ <- function(input_string) {
  gsub("/", "\\\\", input_string, fixed = TRUE)
}

GetImageJParameterPath <- function(ImageJ_app_directory, macro_folder, sub_folder, file_name) {
  file.path(ImageJ_app_directory, "_interface_with_R", "temp", macro_folder, sub_folder, file_name)
}


SplitRawFile <- function(input_file_dir,
                         input_file_name,
                         output_dir_path,
                         ImageJ_app_directory,
                         total_num_slices        = NULL,
                         num_slices_per_substack = 500L,
                         z_step_size             = NULL,
                         splitting_image_text    = "Splitting the following image into .tif files:"
                         ) {

  if (is.null(total_num_slices)) {
    total_num_slices <- 2500L
  }
  if (is.null(z_step_size)) {
    z_step_size <- 3
  }

  input_file_path <- EscapePathForImageJ(file.path(input_file_dir, input_file_name))

  ImageJ_parameter <- paste0(c(output_dir_path, input_file_path, total_num_slices, num_slices_per_substack, z_step_size), collapse = "*")
  ImageJ_parameter <- FormatPathForShell(EscapePathForImageJ(ImageJ_parameter))

  MessageAndLog(paste0(splitting_image_text, " '", input_file_name, "' ..."))
  MessageAndLog(paste0("... which will be placed in the folder: '", output_dir_path, "' ..."))

  old_wd <- setwd(ImageJ_app_directory)

  system_output <- system(paste0("ImageJ-win64.exe --headless -macro Call_from_R__split_for_Terastitcher.ijm", " ", ImageJ_parameter), intern = TRUE)

  setwd(old_wd)
  return(invisible(NULL))
}





# Functions for calling Terastitcher --------------------------------------

terastitcher_command_vec <- c(
  'terastitcher --import --volin="%s" --volin_plugin="TiledXY|3Dseries" --imin_plugin="tiff3D" --ref1=2 --ref2=1 --ref3=3 --vxl1=<<pixel_size>> --vxl2=<<pixel_size>> --vxl3=<<z_step_size>> --projout=xml_import',
  'terastitcher --displcompute --projin="%s"/xml_import.xml --subvoldim=100 --projout=xml_displcomp --sV=50 --sH=50 --sD=25',
  'terastitcher --displproj --projin="%s"/xml_displcomp.xml --projout=xml_displproj',
  'terastitcher --displthres --projin="%s"/xml_displproj.xml --projout=xml_displthres --threshold=0.7',
  'terastitcher --placetiles --projin="%s"/xml_displthres.xml --projout=xml_merging',
  'teraconverter -s="%s/xml_merging.xml" -d="<<full_path>>" --sfmt="TIFF (unstitched, 3D)" --dfmt="TIFF (tiled, 3D)" --resolutions=012'
)



MakeTerastitcherCommands <- function(parent_directory, use_folder_name, pixel_size = 3.26, z_step_size = 3) {
  results_vec <- sprintf(terastitcher_command_vec, use_folder_name)
  results_vec[[length(results_vec)]] <- sub("<<full_path>>", file.path(parent_directory, use_folder_name), results_vec[[length(results_vec)]], fixed = TRUE)
  results_vec[[1]] <- gsub("<<pixel_size>>", as.character(pixel_size), results_vec[[1]], fixed = TRUE)
  results_vec[[1]] <- sub("<<z_step_size>>", as.character(z_step_size), results_vec[[1]], fixed = TRUE)
  return(results_vec)
}





# Functions for logging the progress of Terastitcher  ---------------------

DescribeTerastitcherCommands <- function(desc_vec = terastitcher_command_descriptions) {
  paste0("Terastitcher step ", seq_along(desc_vec), " of ", length(desc_vec), ": ", desc_vec, "...")
}


terastitcher_command_descriptions <- c(
  "Importing parameters",
  "Aligning tiles",
  "Projecting displacements",
  "Tresholding displacements",
  "Placing tiles",
  "Merging tiles"
)


FormatTimeRemaining <- function(time_remaining_line) {
  time_remaining <- sub("TIME REMAINING:\t", "", time_remaining_line)
  time_splits <- strsplit(time_remaining, " minutes and ", fixed = TRUE)[[1]]
  time_splits[[2]] <- sub(" seconds$", "", time_splits[[2]])
  total_seconds <- as.numeric(time_splits[[1]]) * 60 + as.numeric(time_splits[[2]])
  result_string <- paste0(dhms(total_seconds, add_space = TRUE), " remaining...")
  return(result_string)
}

FormatPercentProgress <- function(percent_progress_line) {
  percent_progress <- strsplit(percent_progress_line, "\t", fixed = TRUE)[[1]][[2]]
  percent_numeric <- as.numeric(sub("%", "", percent_progress, fixed = TRUE))
  result_string <- paste0("Progress: ", format(percent_numeric, width = 2), "%")
  return(result_string)
}


TeraStitcherToConsole <- function(line, process) {
  if ((length(current_output) > 1) && (line == "")) {
    print_output <- c(current_output, "")
    if (current_output[[1]] == "") {

      captured_output_list <- c(captured_output_list, list(c(current_output, "")))
      assign("captured_output_list", captured_output_list, envir = globalenv())

      current_output <- gsub("\f", "", current_output, fixed = TRUE)
      is_time_elapsed <- grepl("^Time elapsed: ", current_output)
      if (any(is_time_elapsed)) {
        seconds_elapsed <- sub("Time elapsed: ", "", current_output[is_time_elapsed], fixed = TRUE)
        seconds_elapsed <- sub(" seconds", "", seconds_elapsed, fixed = TRUE)
        seconds_elapsed <- as.numeric(seconds_elapsed)
        if (seconds_elapsed >= 60) {
          print_output <- paste0("Time elapsed: ", dhms(seconds_elapsed, add_space = TRUE))
        } else {
          print_output <- current_output[is_time_elapsed]
        }
        print_output <- c(print_output, "")
      } else if (current_output[[2]] == "OPERATION:\tPairwise displacement computation") {
        if (current_output[[3]] == "PHASE:\t\tInitializing...") {
          print_output <- "Pairwise displacement computation is initializing..."
        } else if (current_output[[3]] == "PHASE:\t\tEnded!") {
          print_output <- "Pairwise displacement computation is complete!"
        } else {
          computation_number_string <- sub("PHASE:\t\tDisplacement computation ", "", current_output[[3]], fixed = TRUE)
          computation_number_splits <- strsplit(computation_number_string, " of ", fixed = TRUE)[[1]]
          print_output <- paste0("Pairwise displacement computation ",
                                 format(as.numeric(computation_number_splits[[1]]), width = nchar(computation_number_splits[[2]])), " of ", computation_number_splits[[2]],
                                 "... ", FormatTimeRemaining(current_output[[4]]), " ", FormatPercentProgress(current_output[[5]])
                                 )
        }
      } else if (grepl("^in VolumeConverter::generateTilesVaa3DRaw", current_output[[2]]) &&
                 (current_output[[3]] == "OPERATION:\tMultiresolution tile generation") &&
                 (current_output[[4]] == "PHASE:\t\tInitializing...")
                 ) {
        print_output <- c(current_output[[2]], "Multiresolution tile generation is initializing...")
      } else if ((current_output[[2]] == "OPERATION:\tMultiresolution tile generation") &&
                 (grepl("^PHASE:\t\t(Generating slices from|Generating resolution|Saving to disc resolution) ", current_output[[3]]))
                 ) {
          current_operation_string <- sub("PHASE:\t\t", "", current_output[[3]], fixed = TRUE)
          current_operation_string <- sub(" og ", " of ", current_operation_string)
          current_operation_string <- paste0(current_operation_string, "... ", paste0(rep(" ", times = 43 - nchar(current_operation_string)), collapse = ""))
          print_output <- paste0("Multiresolution tile generation: ", current_operation_string,
                                 FormatTimeRemaining(current_output[[4]]), " ", FormatPercentProgress(current_output[[5]])
                                 )
      }
    }
    # print_output <- gsub("\f", "", print_output, fixed = TRUE)
    for (my_line in print_output) {
      cat(my_line, fill = TRUE)
    }

    printed_output_list <- c(printed_output_list, list(print_output))
    assign("printed_output_list", printed_output_list, envir = globalenv())
    assign("current_output", "", envir = globalenv())
  } else {
    current_output <- c(current_output, line)
    assign("current_output", current_output, envir = globalenv())
  }
}


DoTerastitch <- function(parent_directory, folder_name, pixel_size = 3.26, z_step_size = 3, allow_sparse_data = FALSE) {
  # old_wd <- setwd(parent_directory)
  commands_vec <- MakeTerastitcherCommands(parent_directory, folder_name, pixel_size = pixel_size, z_step_size = z_step_size)
  if (allow_sparse_data) {
    commands_vec[[1]] <- paste0(commands_vec[[1]], " --sparse_data")
  }
  desc_vec <- DescribeTerastitcherCommands()
  stopifnot(length(commands_vec) == length(desc_vec))
  for (i in seq_along(commands_vec)) {
    MessageAndLog(desc_vec[[i]])
    command_splits <- strsplit(commands_vec[[i]], " ", fixed = TRUE)[[1]]
    command_line_output <- processx:::run(command_splits[[1]],
                                          args = command_splits[2:length(command_splits)],
                                          wd = parent_directory,
                                          windows_verbatim_args = TRUE,
                                          stdout_line_callback = TeraStitcherToConsole,
                                          stderr_line_callback = function(line, process) cat(line)
                                          )
    if (command_line_output[["stderr"]] != "") {
      MessageAndLog(paste0("Terastitcher produced an error: ", sub("\\r\\n$", "", command_line_output[["stderr"]])), send_to_console = FALSE)
    }
    if (grepl("^ERROR: ", command_line_output[["stdout"]])) {
      MessageAndLog(paste0("Terastitcher ", sub("\\r\\n$", "", command_line_output[["stdout"]])), send_to_console = FALSE)
    }
    if (i == 2) {
      LogTime()
    }
  }
  # setwd(old_wd)
  return(invisible(NULL))
}





# Top-level functions -----------------------------------------------------

RunStitch <- function(emission_filters              = NULL,
                      raw_files_root_path           = NULL,
                      split_and_stitch_root_path    = NULL,
                      skip_already_present_splits   = FALSE,
                      skip_already_present_stitches = FALSE,
                      use_metadata                  = TRUE,
                      check_metadata                = TRUE,
                      use_physical_layout_mat       = physical_layout_mat,
                      displacement_in_microns       = 5000,
                      num_slices_per_substack       = 500L,
                      allow_sparse_data             = FALSE,
                      no_rotation                   = FALSE
                      ) {
  if (is.null(raw_files_root_path)) {
    raw_files_root_path <- choose.dir("", caption = "Select a folder or directory of folders of .raw files...")
  }
  raw_files_root_path <- TidyPaths(raw_files_root_path)
  if (is.null(emission_filters)) {
    emission_filter_message <- ""
  } else {
    num_filters <- length(emission_filters)
    if (num_filters > 1) {
      emission_filter_message <- paste0(paste0(emission_filters[seq_len(num_filters - 1)], collapse = ", "), " and ", emission_filters[num_filters])
    } else {
      emission_filter_message <- emission_filters
    }
    emission_filter_message <- paste0("for the ", emission_filter_message, " nm filter", if (num_filters > 1) "s" else "")
  }
  SplitAndTerastitch_args <- list(
    "raw_files_folder_path"         = NA_character_,
    "split_and_stitch_root_path"    = split_and_stitch_root_path,
    "emission_filters"              = emission_filters,
    "skip_already_present_splits"   = skip_already_present_splits,
    "skip_already_present_stitches" = skip_already_present_stitches,
    "use_metadata"                  = use_metadata,
    "check_metadata"                = check_metadata,
    "displacement_in_microns"       = displacement_in_microns,
    "use_physical_layout_mat"       = use_physical_layout_mat,
    "num_slices_per_substack"       = num_slices_per_substack,
    "allow_sparse_data"             = allow_sparse_data,
    "no_rotation"                   = no_rotation
  )
  if (ContainsRawFiles(raw_files_root_path)) {
    MessageAndLog(paste0("Processing the .raw files in the directory: '", raw_files_root_path, "' ..."))
    DisplayPathAsPlot(raw_files_root_path)
    SplitAndTerastitch_args[["raw_files_folder_path"]] <- raw_files_root_path
    result_status <- do.call(SplitAndTerastitch, SplitAndTerastitch_args)
    if (any(result_status)) {
      MessageAndLog(paste0("Stitching of all files in '", raw_files_root_path, "' ", emission_filter_message, " has completed."))
      LogBlankLines(1)
      LogTime()
    } else {
      MessageAndLog(paste0("All files in '", raw_files_root_path, "' were already found to be stitched and split (", emission_filter_message, "), so no processing was done."))
    }
    LogBlankLines(2)
  } else {
    folder_paths_vec <- FindFoldersWithRawFiles(raw_files_root_path)
    if (length(folder_paths_vec) == 0) {
      LogAndStop("No .raw files were found in the chosen folder, nor within any of its subfolders!")
    }
    MessageAndLog(paste0(".raw files in the following folder", if (length(folder_paths_vec) > 1) "s" else "", " will be processed:"))
    for (folder_path in folder_paths_vec) {
      MessageAndLog(paste0("'", folder_path, "'"))
    }
    LogBlankLines(1)
    DisplayPathAsPlot(raw_files_root_path)
    processing_performed <- FALSE
    for (folder_path in folder_paths_vec) {
      MessageAndLog(paste0("Processing the .raw files in the folder: '", folder_path, "' ..."))
      SplitAndTerastitch_args[["raw_files_folder_path"]] <- folder_path
      result_status <- do.call(SplitAndTerastitch, SplitAndTerastitch_args)
      if (any(result_status)) {
        processing_performed <- TRUE
      }
      MessageAndLog(paste0("Processing of the .raw files in the folder '", folder_path, "' has completed."))
      LogBlankLines(2)
    }
    if (processing_performed) {
      MessageAndLog(paste0("Stitching of all subfolders within '", raw_files_root_path, "' ", emission_filter_message, " has completed."))
      LogBlankLines(1)
      LogTime()
    } else {
      MessageAndLog(paste0("All files in '", raw_files_root_path, "' were already found to be stitched and split (", emission_filter_message, "), so no processing was done."))
    }
    LogBlankLines(2)
  }
}



DisplayPathAsPlot <- function(split_and_stitch_root_path) {
  ### The path is displayed in the plots because the output from Terastitcher can crowd out all previous text output on the console,
  ### thus making it difficult to determine which R window is processing which directory (e.g. which external hard disk).
  directory_expression <- as.expression(bquote(bold(.(gsub("/", "\\\\", split_and_stitch_root_path)))))
  plot(1, type = "n", ann = FALSE, axes = FALSE)
  legend("center",
         bty       = "n",
         legend    = sapply(list(expression("Processing and stitching files in this directory:"), directory_expression), function(x) x),
         y.intersp = 1.5,
         cex       = 1.5
         )
  return(invisible(NULL))
}



SplitAndTerastitch <- function(raw_files_folder_path,
                               split_and_stitch_root_path    = NULL,
                               emission_filters              = NULL,
                               skip_already_present_splits   = FALSE,
                               skip_already_present_stitches = FALSE,
                               use_metadata                  = TRUE,
                               check_metadata                = TRUE,
                               displacement_in_microns       = 5000,
                               use_physical_layout_mat       = physical_layout_mat,
                               num_slices_per_substack       = 500L,
                               allow_sparse_data             = TRUE,
                               no_rotation                   = FALSE
                               ) {
  stitch_mat <- PrepareForSplitting(raw_files_folder_path,
                                    split_and_stitch_root_path,
                                    emission_filters = emission_filters,
                                    displacement_in_microns = displacement_in_microns,
                                    use_physical_layout_mat = use_physical_layout_mat,
                                    no_rotation = no_rotation
                                    )
  if (use_metadata) {
    stitch_df <- MetaDataForStitchMat(stitch_mat, check_metadata = check_metadata, tsv_path = MetaDataTSVPath(raw_files_folder_path, split_and_stitch_root_path))
  }
  split_indices <- seq_len(nrow(stitch_mat))
  num_split <- nrow(stitch_mat)
  if (skip_already_present_splits) {
    num_slices_per_substack_vec <- rep(num_slices_per_substack, num_split)
    if (use_metadata) {
      are_already_present_splits <- AlreadyPresentSplits(stitch_mat[, "New_folder_path"],
                                                         total_num_slices_vec_or_list        = stitch_df[, "number_of_z_planes"],
                                                         num_slices_per_substack_vec_or_list = num_slices_per_substack_vec,
                                                         z_step_size_vec_or_list             = stitch_df[, "z_step_size_in_microns"]
                                                         )
    } else {
      are_already_present_splits <- AlreadyPresentSplits(stitch_mat[, "New_folder_path"], num_slices_per_substack_vec_or_list = num_slices_per_substack_vec)
    }
    split_indices <- split_indices[!(are_already_present_splits)]
    if (all(are_already_present_splits)) {
      num_split <- 0
      MessageAndLog(paste0("Split .tif files were already present for all .raw files selected from '", stitch_mat[1, "Old_folder_path"],
                           "' (and were found in '", stitch_mat[1, "New_top_path"], "'), so this step was skipped!"
                           )
                    )
    } else if (any(are_already_present_splits)) {
      num_split <- sum(!(are_already_present_splits))
      MessageAndLog(paste0("Split .tif files were already present for ", sum(are_already_present_splits), " out of ",
                           nrow(stitch_mat), " .raw files selected from '", stitch_mat[1, "Old_folder_path"],
                           "' (and were found in '", stitch_mat[1, "New_top_path"], "'), so these will be skipped!"
                           )
                    )
    }
  }
  for (i in seq_along(split_indices)) {
    index <- split_indices[[i]]
    use_splitting_image_text <- paste0("Splitting image ", i, " of ", length(split_indices), " into .tif files:")
    SplitRawFile(raw_files_folder_path,
                 stitch_mat[index, "Old_file_name"],
                 stitch_mat[index, "New_folder_path"],
                 ImageJ_app_directory,
                 total_num_slices        = if (use_metadata) stitch_df[index, "number_of_z_planes"] else NULL,
                 num_slices_per_substack = num_slices_per_substack,
                 z_step_size             = if (use_metadata) stitch_df[index, "z_step_size_in_microns"] else NULL,
                 splitting_image_text    = use_splitting_image_text
                 )
    if (((i %% 8) == 0) && (i != max(split_indices))) {
      LogTime()
    }
  }
  if (num_split != 0) {
    MessageAndLog(paste0(num_split, " .raw files in '", raw_files_folder_path, "' were split into .tif files."))
    LogTime()
    LogBlankLines(1)
  }
  new_top_path <- stitch_mat[1, "New_top_path"]
  are_unique <- !(duplicated(stitch_mat[, c("Old_folder_path", "New_top_path", "Emission_filter", "Sample_rotation")], MARGIN = 1))
  anything_stitched <- FALSE
  for (i in which(are_unique)) {
    filter <- stitch_mat[i, "Emission_filter"]
    rotation <- stitch_mat[i, "Sample_rotation"]
    filter_folder_path <- file.path(new_top_path, paste0("filter_", filter))
    if (skip_already_present_stitches && AlreadyPresentStitches(file.path(filter_folder_path, paste0("rot_", rotation)))) {
      MessageAndLog(paste0("Stitched files already seemed to be present for the ", filter, " nm emission filter and a sample rotation of ", rotation, "°, so this stitching step was skipped!"))
    } else {
      MessageAndLog(paste0("Calling Terastitcher for the images generated using the ", filter, " nm emission filter and a sample rotation of ", rotation, "°..."))
      if (use_metadata) {
        DoTerastitch(filter_folder_path,
                     paste0("rot_", rotation),
                     pixel_size = stitch_df[i, "pixel_size_in_microns"],
                     z_step_size = stitch_df[i, "z_step_size_in_microns"],
                     allow_sparse_data = allow_sparse_data
                     )
      } else {
        DoTerastitch(filter_folder_path,
                     paste0("rot_", rotation),
                     allow_sparse_data = allow_sparse_data
                     )
      }
      MessageAndLog(paste0("Stitching for the files in '", file.path(filter_folder_path, paste0("rot_", rotation)), "' has completed."))
      anything_stitched <- TRUE
      LogTime()
      LogBlankLines(1)
    }
  }
  result_status <- c("anything_split" = num_split > 0, "anything_stitched" = anything_stitched)
  return(invisible(result_status))
}






# Functions for checking which output files are already present -----------

CheckSplitTIFFsForCompleteness <- function(folder_path,
                                           min_file_size_in_GB        = 0.1,
                                           check_for_known_file_sizes = TRUE,
                                           use_known_file_sizes       = known_file_sizes,
                                           total_num_slices           = NULL,
                                           num_slices_per_substack    = NULL,
                                           z_step_size                = NULL
                                           ) {
  if (is.null(total_num_slices)) {
    total_num_slices <- 2500L
  }
  if (is.null(num_slices_per_substack)) {
    num_slices_per_substack <- 500L
  }
  if (is.null(z_step_size)) {
    z_step_size <- 3
  }
  num_substacks <- ceiling(total_num_slices / num_slices_per_substack)

  file_names <- list.files(folder_path)
  expected_names <- formatC(as.integer(round((seq_len(num_substacks) - 1) * z_step_size * num_slices_per_substack * 10L)), width = 6, flag = "0")
  full_file_names <- paste0(expected_names, ".tif")
  if (!(identical(file_names, full_file_names))) {
    LogAndStop(paste0("The file names in '", folder_path, "' did not conform to the expectation -- i.e., 5 .tif files named: ", paste0(expected_names, collapse = ", "), "!"))
  }
  file_sizes <- file.size(file.path(folder_path, full_file_names)) / 1024 / 1024 / 1024
  rounded_file_sizes <- round(file_sizes, digits = 3)
  are_too_small <- file_sizes < min_file_size_in_GB
  file_name_and_size <- paste0("'", full_file_names, "' (", rounded_file_sizes, " GiB)")
  if (any(are_too_small)) {
    LogAndStop(paste0("The file size for the following files in '", folder_path, "' were smaller than the minimum of ", min_file_size_in_GB, " gigabytes: ",
                      paste0(file_name_and_size[are_too_small], collapse = ", ")
                      )
               )
  }
  if (check_for_known_file_sizes) {
    are_correct_size <- rounded_file_sizes %in% use_known_file_sizes
    if (!(all(are_correct_size))) {
      LogAndStop(paste0("The following files in '", folder_path, "' had an unexpected file size: ",
                        paste0(file_name_and_size[!(are_correct_size)], collapse = ", ")
                        )
                 )
    }
  }
  return(invisible(TRUE))
}


AlreadyPresentSplits <- function(new_folders_vec,
                                 total_num_slices_vec_or_list        = lapply(new_folders_vec, function(x) NULL),
                                 num_slices_per_substack_vec_or_list = lapply(new_folders_vec, function(x) NULL),
                                 z_step_size_vec_or_list             = lapply(new_folders_vec, function(x) NULL)
                                 ) {
  any_files_present <- vapply(new_folders_vec, function(x) length(list.files(x)) != 0, logical(1))
  results_vec <- rep(FALSE, length(new_folders_vec))
  if (any(any_files_present)) {
    are_complete <- vapply(which(any_files_present),
                           function(x) CheckSplitTIFFsForCompleteness(new_folders_vec[[x]],
                                                                      total_num_slices        = total_num_slices_vec_or_list[[x]],
                                                                      num_slices_per_substack = num_slices_per_substack_vec_or_list[[x]],
                                                                      z_step_size             = z_step_size_vec_or_list[[x]]
                                                                      ),
                           logical(1)
                           )
    results_vec[any_files_present] <- are_complete
  }
  return(results_vec)
}


AlreadyPresentStitches <- function(hemisphere_path) {
  present_dirs <- list.dirs(hemisphere_path, full.names = FALSE, recursive = FALSE)
  if (any(grepl("^RES\\(", present_dirs))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}





# Functions for collecting and renaming the stitched files ----------------

CollectHemispheres <- function(stitched_folder_path, stitched_folder_prefix = "R_STITCH_") {
  filter_folder_names <- list.dirs(stitched_folder_path, full.names = FALSE, recursive = FALSE)
  if (length(filter_folder_names) == 0) {
    LogAndStop(paste0("No folders were found in the given directory (specified by the 'stitched_folder_path' argument) in CollectHemispheres: '", stitched_folder_path, "'"))
  }
  hemisphere_mat <- do.call(rbind,
                            lapply(filter_folder_names,
                                   function(x) {
                                     rotation_folder_names <- list.dirs(file.path(stitched_folder_path, x), full.names = FALSE, recursive = FALSE)
                                     results_mat <- cbind("Hemisphere_path" = file.path(stitched_folder_path, x, rotation_folder_names),
                                                          "Filter_folder"   = x,
                                                          "Rotation_folder" = rotation_folder_names
                                                          )

                                   })
                            )
  hemisphere_mat <- cbind(
    hemisphere_mat,
    file_prefix = paste0(sub("^filter_", "em", hemisphere_mat[, "Filter_folder"]),
                           "__",
                           sub(stitched_folder_prefix, "", basename(stitched_folder_path), fixed = TRUE),
                           "__",
                           c("rot_0" = "left", "rot_180" = "right")[hemisphere_mat[, "Rotation_folder"]]
                           )

  )
  return(hemisphere_mat)
}



CollectStitchedFiles <- function(hemisphere_path) {
  folder_names <- list.dirs(hemisphere_path, full.names = FALSE, recursive = FALSE)
  stitched_folders <- grep("^RES\\(", folder_names, value = TRUE)
  if (length(stitched_folders) == 0) {
    MessageAndLog(paste0("No stitched files were found in: '", hemisphere_path, "'!"))
    return(NULL)
  } else {
    resolution_splits <- strsplit(sub("\\)", "", sub("^RES\\(", "", stitched_folders)), "x", fixed = TRUE)
    tiff_paths <- vapply(stitched_folders, function(x) {
      sub_folder <- list.dirs(file.path(hemisphere_path, x), full.names = FALSE, recursive = FALSE)
      sub_sub_folder <- list.dirs(file.path(hemisphere_path, x, sub_folder), full.names = FALSE, recursive = FALSE)
      sub_sub_path <- file.path(hemisphere_path, x, sub_folder, sub_sub_folder)
      tiff_file <- list.files(file.path(sub_sub_path))
      if (!((length(tiff_file) == 1) && (grepl("\\.tif$", tiff_file)))) {
        LogAndStop(paste0("The folder '", sub_sub_path, "' did not contain the expected 1 .tif file!"))
      }
      return(file.path(sub_sub_path, tiff_file))
    }, "")
    results_df <- data.frame(
      "Res_number" = NA_integer_,
      "tif_path"   = tiff_paths,
      "x_res"      = as.integer(sapply(resolution_splits, "[", 2)),
      "y_res"      = as.integer(sapply(resolution_splits, "[", 1)),
      "z_res"      = as.integer(sapply(resolution_splits, "[", 3)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    results_df[, "Res_number"] <- order(results_df[, "x_res"], results_df[, "y_res"], results_df[, "x_res"])
    results_df <- results_df[results_df[, "Res_number"], ]
    rownames(results_df) <- NULL
    return(results_df)
  }
}




PrepareForCopying <- function(stitched_folder_path, stitched_folder_prefix = "R_STITCH_") {
  assign("delete_stitched_folder_path", stitched_folder_path, envir = globalenv())
  hemisphere_mat <- CollectHemispheres(stitched_folder_path, stitched_folder_prefix = stitched_folder_prefix)
  hemisphere_df_list <- lapply(seq_len(nrow(hemisphere_mat)), function(x) {
    res_df <- CollectStitchedFiles(hemisphere_mat[x, "Hemisphere_path"])
    if (is.null(res_df)) {
      return(NULL)
    } else {
      results_df <- data.frame(res_df,
                               hemisphere_mat[rep(x, nrow(res_df)), ],
                               stringsAsFactors = FALSE
                               )
      return(results_df)
    }
  })
  hemisphere_mat_list <- Filter(Negate(is.null), hemisphere_df_list)
  hemisphere_df <- do.call(rbind.data.frame, c(hemisphere_df_list, list(stringsAsFactors = FALSE)))

  filter_vec <- sub("^filter_", "", hemisphere_df[, "Filter_folder"])
  
  num_res <- length(unique(hemisphere_df[, "Res_number"]))
  subfolder_names_vec <- c("2) Low", "1) Mid", "0) Full")
  file_res_names_vec <- paste0(c("Lo", "Mid", "Full"), "Res")
  if (num_res == 2) {
    subfolder_names_vec <- subfolder_names_vec[-2]
    file_res_names_vec <- file_res_names_vec[-2]
  } else if (num_res != 3) {
    stop("The 'PrepareForCopying' function only supports 2 or 3 different resolutions!")
  }
  
  hemisphere_df[, "New_subfolder_name"] <- paste0(subfolder_names_vec[hemisphere_df[, "Res_number"]],
                                                  "_resolution_",
                                                  filter_vec,
                                                  "_filter" # _nm_filter
                                                  )
  hemisphere_df[, "Full_new_name"] <- paste0(file_res_names_vec[hemisphere_df[, "Res_number"]],
                                             "_",
                                             hemisphere_df[, "file_prefix"],
                                             "_x",
                                             hemisphere_df[, "x_res"],
                                             "_y",
                                             hemisphere_df[, "y_res"],
                                             "_z",
                                             hemisphere_df[, "z_res"],
                                             ".tif"
                                             )
  return(hemisphere_df)
}


GetRootPath <- function(selected_path) {
  if (ContainsRawFiles(selected_path)) {
    root_path <- dirname(selected_path)
  } else {
    root_path <- selected_path
  }
  root_path <- TidyPaths(root_path)
  return(root_path)
}


CompileAllStitchedFiles <- function(selected_path,
                                    output_folder_path         = NULL,
                                    old_stitched_folder_prefix = "R_STITCH_",
                                    rename_stitched_folder     = FALSE,
                                    new_stitched_folder_prefix = "STITCHED_",
                                    skip_present_files         = FALSE,
                                    make_sub_folders           = TRUE
                                    ) {

  root_path <- GetRootPath(selected_path)

  if (is.null(output_folder_path)) {
    output_folder_path <- file.path(root_path, "Stitched")
  }
  output_folder_path <- TidyPaths(output_folder_path)

  top_folder_names <- list.dirs(root_path, full.names = FALSE, recursive = FALSE)
  are_stitched_folders <- grepl(paste0("^", old_stitched_folder_prefix), top_folder_names)
  if (!(any(are_stitched_folders))) {
    LogAndStop(paste0("No folders were found in the selected directory that met the criteria for the folder name (i.e. beginning with '", old_stitched_folder_prefix, "')!"))
  }
  if (!(dir.exists(output_folder_path))) {
    MessageAndLog(paste0("The output folder '", output_folder_path, "' did not exist and was created!"))
    dir.create(output_folder_path)
  }
  stitched_folder_paths <- file.path(root_path, top_folder_names[are_stitched_folders])
  stitched_folder_paths <- TidyPaths(stitched_folder_paths)
  any_copied <- FALSE
  for (stitched_folder_path in stitched_folder_paths) {
    brain_df <- PrepareForCopying(stitched_folder_path, stitched_folder_prefix = old_stitched_folder_prefix)
    if (make_sub_folders) {
      for (new_sub_folder in unique(brain_df[, "New_subfolder_name"])) {
        dir.create(TidyPaths(file.path(output_folder_path, new_sub_folder)), showWarnings = FALSE)
      }
    }
    if (make_sub_folders) {
      destination_vec <- file.path(output_folder_path, brain_df[, "New_subfolder_name"], brain_df[, "Full_new_name"])
    } else {
      destination_vec <- file.path(output_folder_path, brain_df[, "Full_new_name"])
    }
    destination_vec <- TidyPaths(destination_vec)
    destination_exists_vec <- file.exists(destination_vec)
    if (skip_present_files && all(destination_exists_vec)) {
      MessageAndLog(paste0("All stitched files in '", stitched_folder_path, "' were already found in '", output_folder_path, "', so these files were not copied."))
    } else {
      num_to_copy <- sum(!(destination_exists_vec))
      if (skip_present_files && any(destination_exists_vec)) {
        num_exist <- sum(destination_exists_vec)
        MessageAndLog(paste0("The following", if (num_exist > 1) paste0(num_exist, " ") else "", " stitched file",
                             if (num_exist > 1) "s were" else " was", " already present in '", output_folder_path,
                             "', so these will be skipped: "
                             )
                      )
        for (i in which(destination_exists_vec)) {
          MessageAndLog(paste0("- '", brain_df[i, "Full_new_name"], "'"))
        }
      }
      MessageAndLog(paste0("Copying ", num_to_copy, " stitched file", if (num_to_copy > 1) "s" else "", " from '", stitched_folder_path, "' to '", output_folder_path, "' ..."))
      for (i in seq_len(nrow(brain_df))) {
        if (skip_present_files && destination_exists_vec[[i]]) {
          next
        } else {
          if (destination_exists_vec[[i]]) {
            MessageAndLog(paste0("A file already existed at '", destination_vec[[i]], "' and was deleted."))
            unlink(destination_vec[[i]])
          }
          MessageAndLog(paste0("Copying the file with the new name: '", brain_df[i, "Full_new_name"], "' ..."))
          shell(paste0(c("copy", FormatPathForShell(brain_df[i, "tif_path"]), FormatPathForShell(destination_vec[[i]])), collapse = " ")) # Use a system call for copying the file, because file.copy is far too slow...
          any_copied <- TRUE
        }
      }
      new_stitched_folder_name <- sub(paste0("^", old_stitched_folder_prefix), new_stitched_folder_prefix, basename(stitched_folder_path))
      if (rename_stitched_folder) {
        file.rename(from = stitched_folder_path, to = TidyPaths(file.path(dirname(stitched_folder_path), new_stitched_folder_name)))
      }
      LogTime()
      LogBlankLines(1)
    }
  }
  if (any_copied) {
   MessageAndLog(paste0("Stitched files have been copied from the '", root_path, "' directory to '", output_folder_path, "'."))
  }
}



CompileStitchedFromTo <- function(from_path = NULL, to_path = NULL, stitched_folder_prefix = "R_STITCH_", ...) {
  if (is.null(from_path)) {
    from_path <- choose.dir("", caption = paste0("Take files from... (must contain folders beginning with '", stitched_folder_prefix, "')"))
  }
  if (is.null(to_path)) {
    to_path <- choose.dir("", caption = "Copy to...")
  }
  CompileAllStitchedFiles(from_path, output_folder_path = to_path, old_stitched_folder_prefix = stitched_folder_prefix, ...)
  return(invisible(NULL))
}




# Functions for updating a log file ---------------------------------------

SetUpLogFile <- function(root_directory, log_sub_folder = "Stitch_log") {
  if (!(is.null(log_file_path))) {
    stop(paste0("The global variable 'log_file_path' already existed, with the value: '", log_file_path, "'!"))
  }
  log_file_name <- paste0("StitchLog_", format(Sys.time(), format = "%Y-%m-%d_%H%M"), ".txt")
  if (!(is.null(log_sub_folder))) {
    dir.create(TidyPaths(file.path(root_directory, log_sub_folder)), showWarnings = FALSE)
    log_file_path <- file.path(root_directory, log_sub_folder, log_file_name)
  } else {
    log_file_path <- file.path(root_directory, log_file_name)
  }
  log_file_path <- TidyPaths(log_file_path)
  assign("log_file_path", log_file_path, envir = globalenv())
  initial_message <- paste0("This log file was generated by Lukas' script for stitching .raw files, and was created on ",
                            format(Sys.time(), format = "%Y-%m-%d at %H:%M:%S.")
                            )
  initial_message <- c(initial_message, "", "")
  initial_message <- as.matrix(initial_message)
  log_time_zero <- Sys.time()
  assign("log_time_zero", log_time_zero, envir = globalenv())
  write.table(initial_message, file = log_file_path, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  return(invisible(NULL))
}



MessageAndLog <- function(message_string, create_log_file = TRUE, stop = FALSE, send_to_console = TRUE) {
  ### Requires the global variable "log_file_path" ###
  if (!(stop) && send_to_console) {
    message(message_string)
  }
  if (create_log_file && !(is.null(log_file_path))) {
    old_log_contents <- as.matrix(read.table(file.path(log_file_path), quote = "", header = FALSE, row.names = NULL,
                                             sep = "\t", stringsAsFactors = FALSE, blank.lines.skip = FALSE
                                             )
                                  )
    new_log_contents <- rbind(old_log_contents, message_string)
    write.table(new_log_contents, file = log_file_path, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t", append = FALSE)
  }
  if (stop) {
    stop(message_string)
  }
  return(invisible(NULL))
}


LogAndStop <- function(stop_string, create_log_file = TRUE) {
  MessageAndLog(message_string = stop_string, create_log_file = create_log_file, stop = TRUE)
}

LogBlankLines <- function(num_lines = 1L, create_log_file = TRUE) {
  if (create_log_file && !(is.null(log_file_path))) {
    for (i in seq_len(num_lines)) {
      MessageAndLog("", create_log_file = create_log_file)
    }
  }
}


LogTime <- function() {
  if (!(is.null(log_time_zero) || is.null(log_file_path))) {
    time_difference <- as.double(Sys.time() - log_time_zero, units = "secs")
    time_string <- dhm(time_difference)
    message_string <- paste0(time_string, " have elapsed since logging began.")
    MessageAndLog(message_string)
  }
}





# Example commands (commented out) ----------------------------------------

# CompileAllStitchedFiles(GetRootPath(raw_files_folder), skip_present_files = TRUE, rename_stitched_folder = FALSE)
# MetaDataForStitchMat(PrepareForSplitting(raw_files_folder))







# Process images ----------------------------------------------------------







# ### Prepare for stitching ###
#
# stitched_and_renamed_folder <- split_and_stitch_folder
#
# SetUpLogFile(GetRootPath(split_and_stitch_folder))
#
# CheckAllForCompleteness(raw_files_folder, split_and_stitch_folder, raise_metadata_error = FALSE, raise_error = FALSE, emission_filters = c("405D", "594"),
#                         no_rotation = TRUE, use_physical_layout_mat = cbind(2:1)
#                         )
#
#
#
# ### Perform stitching ###
#
# RunStitch(raw_files_root_path = raw_files_folder, split_and_stitch_root_path = split_and_stitch_folder, emission_filters = "405D",
#           skip_already_present_splits = TRUE, skip_already_present_stitches = TRUE, check_metadata = FALSE, allow_sparse_data = FALSE,
#           no_rotation = TRUE, use_physical_layout_mat = cbind(2:1)
#           )
# RunStitch(raw_files_root_path = raw_files_folder, split_and_stitch_root_path = split_and_stitch_folder, emission_filters = "594",
#           skip_already_present_splits = TRUE, skip_already_present_stitches = TRUE, check_metadata = FALSE, allow_sparse_data = FALSE,
#           no_rotation = TRUE, use_physical_layout_mat = cbind(2:1)
#           )
#
#
#
# ### Copy all stitched files (with new, sensible names) into one folder with the name 'Stitched'
#
# CompileStitchedFromTo(split_and_stitch_folder, rename_and_copy_folder, skip_present_files = TRUE)






# Wren study --------------------------------------------------------------


### Get user input on the input and output directories ###

raw_files_folder        <- choose.dir("", caption = "Select a folder or directory of folders of .raw files...")
split_and_stitch_folder <- choose.dir("", caption = "Select a folder for stitching...")
rename_and_copy_folder  <- choose.dir("", caption = "Select a folder to copy stitched files...")



### Prepare for stitching ###

stitched_and_renamed_folder <- split_and_stitch_folder

SetUpLogFile(GetRootPath(split_and_stitch_folder))

# CheckAllForCompleteness(raw_files_folder, split_and_stitch_folder, raise_metadata_error = TRUE, raise_error = TRUE, emission_filters = c("498", "565"))



### Perform stitching ###

RunStitch(raw_files_root_path = raw_files_folder, split_and_stitch_root_path = split_and_stitch_folder, emission_filters = "515LP",
          skip_already_present_splits = TRUE, skip_already_present_stitches = TRUE, check_metadata = FALSE, allow_sparse_data = FALSE
          )



### Copy all stitched files (with new, sensible names) into one folder with the name 'Stitched'

CompileStitchedFromTo(split_and_stitch_folder, rename_and_copy_folder, skip_present_files = TRUE)














