#!/bin/bash
source ./env-jdk17.sh
java -version

cd target/classes
java test.FindExtremaAccuracyTest

