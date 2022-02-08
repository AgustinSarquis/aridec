context("Template Structure")

test_that("all entries from the database can be read", {
  expect_silent(loadEntries(path=path_to_data))
})

database<-loadEntries(path=path_to_data)

test_that("entry names correspond to metadata template", {
  template=yaml::yaml.load_file(input=paste(path_to_data, "metadata Template.yaml", sep=""))
  for(i in 1:length(database)){
    entry=database[[i]]
    expect_equal(names(entry)[-(12:13)],names(template))
  }

  for(i in 1:length(database)){
    entry=database[[i]]
    expect_equal(names(entry$siteInfo),names(template$siteInfo))
  }

  # allows new fields to be added, but none to be removed
  for(i in 1:length(database)){
    entry=database[[i]]
    if(length(match(names(template$experimentInfo),names(entry$experimentInfo)))!=length(names(template$experimentInfo)))
      {cat(names(database[i]),"\n")}
    expect_equal(length(match(names(template$experimentInfo),names(entry$experimentInfo))),
                 length(names(template$experimentInfo)))
  }
})

test_that("check first column in timeSeries file is called 'Time' ",{
  for(i in 1:length(database)){
    entry=database[[i]]
    expect_equal(colnames(entry$timeSeries)[1], "Time")

  }
})

test_that("check that timeSeries time variable matches allowable units",{
  for(i in 1:length(database)){
    v1.unit<-database[[i]][["variables"]][[1]]["units"]
   if(length(which(v1.unit == "days" | v1.unit == "weeks"| v1.unit == "months"| v1.unit == "years"))==0) {
     cat(names(database)[i],"\n")
    }
    expect_true(v1.unit == "days" | v1.unit == "weeks"| v1.unit == "months"| v1.unit == "years")
  }
})

test_that("variable names in timeSeries correspond to names in metadata file",{
  for(i in 1:length(database)){
    entry=database[[i]]
    tsName=as.character(colnames(entry$timeSeries))
    varname=as.character(sapply(entry$variables, function(x){x$name}))
    expect_equal(tsName, varname)
  }
})

testthat::test_that("check that each table in variables list has the same number of fields",{
  n.fields <- lapply(database, function(x) {
    len <- length(unique(unlist(lapply(x$variables[-1], function(y) length(y)))))
    return(len)
  })
  for(i in seq_along(n.fields)){
    if(n.fields[[i]]!=1) {cat(names(n.fields[i]),"\n")}
    expect_equal(n.fields[[i]], 1)
  }
})
