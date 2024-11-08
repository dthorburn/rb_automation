IPR_DIR=/mnt/sda1/interproscan

export PATH=/home/resurrect/conda/envs/iprscan/bin:$PATH

cd $IPR_DIR


if [ ! -f IPR_READY ]
then
	wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.71-102.0/interproscan-5.71-102.0-64-bit.tar.gz
	wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.71-102.0/interproscan-5.71-102.0-64-bit.tar.gz.md5
	md5sum -c interproscan-5.71-102.0-64-bit.tar.gz.md5 && touch IPR_READY
fi

if [ -f IPR_READY ]
then
	tar -pxvzf interproscan-5.71-102.0-*-bit.tar.gz
	cd ${IPR_DIR}/interproscan-5.71-102.0
	python3 setup.py -f interproscan.properties
else
	echo "WARN: interproscan.tar.gz checksums failed."
	exit 1
fi
