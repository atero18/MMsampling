make_exclusions_lintr <- function()
{
  # Files and folders completly ignores
  ignored <- file.path(c("inst/exclusions_lintr.R",
                         "vignettes/",
                         "tests/testthat.R",
                         "tests/testthat/helper.R",
                         "R/AHR_single2.R"))

  exclusions_list <- as.list(ignored)

  testFilenames <- list.files(file.path("tests/testthat"), pattern = "^test_")

  # Indentation linter is ignore within each test file
  for (testFile in testFilenames)
  {
    pathFile <- file.path("tests/testthat/", testFile)
    exclusions_list[[pathFile]] <- list(indentation_linter = Inf)
  }

  exampleFilenames <- list.files(file.path("inst/examples"), pattern = ".R$")


  # Remove some linters for example files
  for (exampleFile in exampleFilenames)
  {
    pathFile <- file.path("inst/examples/", exampleFile)
    exclusions_list[[pathFile]] <- list(implicit_integer_linter = Inf,
                                        commented_code_linter = Inf,
                                        undesirable_function_linter = Inf)
  }



  exclusions_list
}
