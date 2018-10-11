#!/bin/bash

curl -F file=@${1} -F output=text -F offcheck=u http://sbml.org/validator/
