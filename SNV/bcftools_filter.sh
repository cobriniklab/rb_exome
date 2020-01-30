#!/bin/bash

vcf=$1

bcftools filter -i 'FILTER="."' $vcf