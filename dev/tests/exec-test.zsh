#!/usr/bin/env zsh

set -e
function check {
}
cd ${1:-.}

source test.t
_TEST_INPUT_DIR=$PWD
_TESTING_DIR=$PWD/testing_dir

function cleanup {
    rm -rf ${_TESTING_DIR}
}
mkdir -p ${_TESTING_DIR}
trap cleanup EXIT INT TERM

cd ${_TESTING_DIR}
for i in ${INPUTS} ${DEPENDENCIES}; cp ${_TEST_INPUT_DIR}/${i} .

_STDOUT_FILE=test_framework.stdout
execute >_STDOUT_FILE
if [ -n "${STDOUT+x}" ]; then
    if ! diff -u ${_TEST_INPUT_DIR}/${STDOUT} _STDOUT_FILE >/dev/null; then
        echo "Stdout was not as expected" >&2
        echo >&2
        diff -u ${_TEST_INPUT_DIR}/${STDOUT} _STDOUT_FILE >&2
        exit 1
    fi
fi

for i in ${INPUTS}; diff -u ${_TEST_INPUT_DIR}/${i} $(basename ${i})
check
for i in ${OUTPUTS}; diff -u ${_TEST_INPUT_DIR}/${i} $(basename ${i})
echo "OK"
