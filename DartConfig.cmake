# Dashboard is opened for submissions for a 24 hour period starting at
# the specified NIGHLY_START_TIME. Time is specified in 24 hour format.
SET (NIGHTLY_START_TIME "23:00:00 PST")

#
# Dart server to submit results (used by client)
#
SET (PROJECT_UNIX_NAME "simbody")
SET (DROP_SITE "simtk.org")
SET (DROP_LOCATION "/dart/${PROJECT_UNIX_NAME}")
SET (DROP_SITE_USER "anonymous")
SET (DROP_SITE_PASSWORD "")
SET (TRIGGER_SITE "http://${DROP_SITE}/dart/cgi-bin/TriggerDart.cgi?project=${PROJECT_UNIX_NAME}")

# Project Home Page
SET (PROJECT_URL "https://${DROP_SITE}/home/${PROJECT_UNIX_NAME}")

#
# Dart server configuration 
#
SET (CVS_WEB_URL "https://${DROP_SITE}/websvn/wsvn/${PROJECT_UNIX_NAME}")

