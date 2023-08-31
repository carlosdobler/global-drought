
library(tidyverse)
library(stars)


https://nex-gddp-cmip6.s3.us-west-2.amazonaws.com/index.html#NEX-GDDP-CMIP6/ACCESS-CM2/historical/r1i1p1f1/hurs/


library(httr)

library(httr)
library(xml2)
url <- "https://nex-gddp-cmip6.s3.us-west-2.amazonaws.com"
url <- "https://nex-gddp-cmip6.s3-us-west-2.amazonaws.com/NEX-GDDP-CMIP6/ACCESS-CM2/historical/"
response <- GET(url)
content <- content(response, as = "text")
parsed_content <- read_xml(content)
parsed_content <- xml_ns_strip(parsed_content)

pattern <- "ACCESS-CM2/historical"

filenames <- xml_text(xml_find_all(parsed_content, "//Contents/Key/"))
filenames <- xml_text(xml_find_all(parsed_content, paste0("//Contents/Key[contains(., '", pattern, "')]")))


marker <- tail(filenames, 1)
url_with_marker <- modify_url(url, query = list(marker = marker))
response <- GET(url_with_marker)
new_content <- content(response, as = "text")
new_parsed_content <- read_xml(new_content)
new_parsed_content <- xml_ns_strip(new_parsed_content)
xml_text(xml_find_all(new_parsed_content, "//Contents/Key"))


contents_nodes <- xml_find_all(parsed_content, "//Contents")

key_nodes <- xml_find_all(parsed_content, "//Key")
filenames <- xml_text(key_nodes)



parsed_content %>% 
  xml_find_chr("//Key")

filenames <- xml_text(xml_find_all(parsed_content, "//Key"))


base_url <- "https://nex-gddp-cmip6.s3-us-west-2.amazonaws.com/NEX-GDDP-CMIP6/ACCESS-CM2/historical/r1i1p1f1/hurs/"
output_dir <- "/path/to/output/directory/"

# Send an HTTP GET request to retrieve the directory listing
response <- GET(base_url)

# Extract the file names from the directory listing
file_names <- gsub("(?<=href=\")(.*?)(?=\")", "", content(response, as = "text"), perl = TRUE)

# Loop through the file names and download each file
for (file_name in file_names) {
  # Construct the file URL
  file_url <- paste0(base_url, file_name)
  
  # Download the file
  GET(file_url, write_disk(paste0(output_dir, file_name)))
}
