## Utilities

### Create Accessory files

The file, `create_files`, is an R script meant to run from the command-line. It accesses the Google Sheet here:

https://docs.google.com/spreadsheets/d/1n3hXzzhrHZgClLD8P3cyIrK_6YgdjtdGlLswNAuoKSI/edit?usp=sharing

And creates a tidy csv file of software after adding a few automatically-calculated columns. It also writes a json-format file of the same data, but with one record per software package. Categories are collapsed to a json array. To run:

```
create_files.R /path/to/git/repo
# or
create_files.R /path/to/directory/for/accessory_files
```
