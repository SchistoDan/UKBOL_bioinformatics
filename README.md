# UKBOL_bioinformatics
SOP for processing raw genome skims generated from museum specimens, recovering barcoding gene(s), and mitogenome assembly.



## 1. Generating sample metadata
This section describes the process of generating the sample_metadata.csv file from specimen metadata.
1. The sample_metadata sheet in the specimen sampling [`UK*_processing.xlsx`](https://naturalhistorymuseum.sharepoint.com/:f:/r/sites/UKBarcodeofLife/Shared%20Documents/UKBOL_accelerated/UKBOL_Collection/batches?csf=1&web=1&e=ciVFcy) will be automatically populated when filling out the 'Processing' sheet of the same `UK*_processing.xlsx` file.
2. Once the speciment sampling spreadsheet is completed, ensure you are 'on' the sample_metadata sheet, click `File`, `Export`, and `Download CSV as UTF-8` (encoding). This will automatically save the sample_metadata sheet to your `Downloads/` folder as `UK*_processing(sample_metadata).csv`.
3. Copy this newly generated `sample_metadata.csv` into [`UKBOL_Bioinf/sample_metadata/`](https://naturalhistorymuseum.sharepoint.com/:f:/r/sites/UKBarcodeofLife/Shared%20Documents/UKBOL_accelerated/UKBOL_Bioinf/sample_metadata?csf=1&web=1&e=ZIbj7G) and update the [`batch tracking`](https://naturalhistorymuseum.sharepoint.com/:x:/r/sites/UKBarcodeofLife/Shared%20Documents/UKBOL_accelerated/UKBOL%20batch%20tracker.xlsx?d=w5ccef39617a3412b86db93ecedacffbe&csf=1&web=1&e=SqdTxL) spreadsheet accordingly.
- The `sample_metadata.csv` will contain the following fields:
| Sample ID | phylum | class | order | family | subfamily | genus | species | subspecies | type_status | specimen_voucher | collection_date | geographic_location | latitude | longitude |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| UK001-A01 | Arthropoda | Insecta | Coleoptera | Carabidae |  | Carabus | speciousus | subspeciousis | Type | NHMUK012345678 |	United Kingdom | |

4. Copy `sample_metadata.csv` from `UKBOL_Bioinf/sample_metadata/` over to the Crop Diversity cluster (`~/projects/nhm/museomix/UKBOL_accelerated/sample_metadata/`) and perform a sanitisation check, correcting any data if necessary. The `sample_metadata.csv` is now ready for downstream analyses.
*Sanitisation check: Manually check ID and taxonomic lineage columns for specieal characters, trailing spaces, tabs, etc and remove them.




## 2. Sample processing
1..
2..
3..
4..
5..

Transfer sequence data from NHM cluster to Crop Diversity cluster:
```
```

`~/projects/nhm/museomix/UKBOL_accelerated/sample_sheets/`

```
```



## 3. Recovering and validating barcodes with [BeeGees](https://github.com/bge-barcoding/BeeGees)

Requirements at a glance:
- 



## 4. Assemblying and validating mitochondrial genomes with (modified) [skim2mito](https://github.com/SchistoDan/skim2mito/tree/main)

Requirements at a glance:
- 






