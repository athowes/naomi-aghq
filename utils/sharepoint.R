#' Setting details required by spud to reach the sharepoint as environmental variables
pass <- "ENTER_PASS_HERE"

Sys.setenv(
  "SHAREPOINT_USERNAME" = "ath19@ic.ac.uk",
  "SHAREPOINT_PASS" = pass,
  "SHAREPOINT_URL" = "https://imperiallondon.sharepoint.com/",
  "SHAREPOINT_SITE" = "HIVInferenceGroup-WP"
)
