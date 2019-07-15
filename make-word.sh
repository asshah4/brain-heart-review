#!/bin/bash
pandoc review.md -s --reference-doc=word_style.docx --bibliography=hrv.bib -o review.docx
