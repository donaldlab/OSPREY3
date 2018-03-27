
# download gdrive tool
if test -s gdrive
then
    echo "Already have gdrive"
else
    wget -O gdrive "https://docs.google.com/uc?id=0B3X9GlR6EmbnQ0FtZmJJUXEyRTA&export=download"
    chmod +x gdrive
fi

# get the build file (just the python build for now)
DIST=$(find ./build/distributions -maxdepth 1 -name "osprey-python-*.zip")
echo "Found distribution file: $DIST"

# add the build number to the filename
FILENAME=$(basename $DIST .zip)-b$TRAVIS_BUILD_NUMBER.zip
echo "Uploading as $FILENAME ..."

./gdrive --refresh-token $GDRIVE_REFRESH_TOKEN upload --parent $GDRIVE_DIR --name $FILENAME "$DIST"

