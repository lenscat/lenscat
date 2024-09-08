#!/usr/bin/bash

# Update output_catalog
echo "Updating output_catalog"
# Run python to get the output
ipython -c 'import lenscat; from astropy.table.pprint import conf; conf.max_width=-1; lenscat.catalog' > output_catalog.out 2>/dev/null
sed -i '1d' output_catalog.out
# Update README.md with the recorded output
awk '/```/{f=0} !f; sub(/```output_catalog/," "){print; f=1}' README.md > README.md.tmp
sed -i -e '/```output_catalog/r output_catalog.out' README.md.tmp
# Replace README.md with the updated README.md.tmp
mv README.md.tmp README.md

# Update output_search_confident_cluster
echo "Updating output_search_confident_cluster"
# Run python to get the output
ipython -c 'import lenscat, numpy; from astropy.table.pprint import conf; conf.max_width=-1; lenscat.catalog.search(grading="confident", lens_type="cluster", zlens_range=(1,numpy.inf)).show_ref()' > output_search_confident_cluster.out 2>/dev/null
sed -i '1d' output_search_confident_cluster.out
# Update README.md with the recorded output
awk '/```/{f=0} !f; sub(/```output_search_confident_cluster/," "){print; f=1}' README.md > README.md.tmp
sed -i -e '/```output_search_confident_cluster/r output_search_confident_cluster.out' README.md.tmp
# Replace README.md with the updated README.md.tmp
mv README.md.tmp README.md

# Update output_search_GW170814
echo "Updating output_search_GW170814"
# Download the data
wget https://dcc.ligo.org/public/0157/P1800381/007/GW170814_skymap.fits.gz
# Run python to get the output
ipython -c 'import lenscat; from astropy.table.pprint import conf; conf.max_width=-1; lenscat.catalog.search(lens_type="galaxy").crossmatch("GW170814_skymap.fits.gz")' > output_search_GW170814.out 2>/dev/null
sed -i '1d' output_search_GW170814.out
# Update README.md with the recorded output
awk '/```/{f=0} !f; sub(/```output_search_GW170814/," "){print; f=1}' README.md > README.md.tmp
sed -i -e '/```output_search_GW170814/r output_search_GW170814.out' README.md.tmp
# Replace README.md with the updated README.md.tmp
mv README.md.tmp README.md

# Update output_search_GRB
echo "Updating output_search_GRB"
# Download the data
wget https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2024/bn240229588/quicklook/glg_healpix_all_bn240229588.fit
# Run python to get the output
ipython -c 'import lenscat; from astropy.table.pprint import conf; conf.max_width=-1; lenscat.catalog.crossmatch("glg_healpix_all_bn240229588.fit")' > output_search_GRB.out 2>/dev/null
sed -i '1d' output_search_GRB.out
# Update README.md with the recorded output
awk '/```/{f=0} !f; sub(/```output_search_GRB/," "){print; f=1}' README.md > README.md.tmp
sed -i -e '/```output_search_GRB/r output_search_GRB.out' README.md.tmp
# Replace README.md with the updated README.md.tmp
mv README.md.tmp README.md

echo "Cleaning up"
# Clean up
rm output_catalog.out
rm output_search_confident_cluster.out
rm GW170814_skymap.fits.gz
rm output_search_GW170814.out
rm glg_healpix_all_bn240229588.fit
rm output_search_GRB.out
