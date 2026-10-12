# Export as HTML

# Create a new environment to load the help objects into
DiffBind_help <- new.env()
lazyLoad(
  file.path(system.file("help", package = "DiffBind"), "DiffBind"),
  envir = DiffBind_help
)

# Function to export help files to an HTML document
export_help_to_html <- function(env, file_path) {
  con <- file(file_path, "w")
  on.exit(close(con))
  
  for (obj_name in ls(env)) {
    help_content <- get(obj_name, envir = env)
    if (inherits(help_content, "Rd")) {
      html_content <- capture.output(tools::Rd2HTML(help_content))
      writeLines(paste("Help for:", obj_name), con)
      writeLines(html_content, con)
    }
    writeLines("\n\n", con)
  }
}

# Export the help files to an HTML document
export_help_to_html(DiffBind_help, "/vast/projects/Immune_Aging/Itokawa_et_al_2022_ATACseq/DiffBind_help.html")
