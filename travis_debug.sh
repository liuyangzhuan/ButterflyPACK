#! /usr/bin/env bash

# <token> can be found at https://travis-ci.com/profile
# <jobid> can be found by expanding "Build system information" in the build log
# after running this script (without sudo), got to build log to find SSH command
# after SSH connection, do
# travis_run_before_install
# travis_run_install
# travis_run_before_script
# travis_run_script
# travis_run_after_success
# travis_run_after_failure
# travis_run_after_script


curl -s -X POST \
-H "Content-Type: application/json" \
-H "Accept: application/json" \
-H "Travis-API-Version: 3" \
-H "Authorization: token SQf8BFJHgY4im2WoT-s-pg" \
-d '{ "quiet": true }' \
https://api.travis-ci.com/job/225747268/debug
