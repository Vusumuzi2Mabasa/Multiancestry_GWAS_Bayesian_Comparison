# Render the shared research narrative; regenerate figures separately if needed.
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)])
script <- gsub('~+~', ' ', script, fixed=TRUE)
root <- normalizePath(file.path(dirname(script), '..'), mustWork=TRUE)
python <- Sys.getenv('PYTHON', unset='python')
status <- system2(python, shQuote(file.path(root, 'scripts', 'build_word_report.py')))
if(status != 0L) stop('Word build failed. Install python-docx and set PANDOC to your Pandoc executable.')
