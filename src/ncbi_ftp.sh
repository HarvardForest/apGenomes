#!/bin/bash

### After configuring and saving the .sh script, make it executable:
### chmod +x ftpscript.sh

HOST=lftp-private.ncbi.nlm.nih.gov  #This is the FTP servers host or IP address.
USER=subftp          #This is the FTP user that has access to the server.
PASS=w4pYB9VQ          #This is the password for the FTP user.

# Call 1. Uses the ftp command with the -inv switches.
#-i turns off interactive prompting.
#-n Restrains FTP from attempting the auto-login feature.
#-v enables verbose and progress.

ftp -inv $HOST <<EOF
EOF

# Call 2. Here the login credentials are supplied by calling the variables.

user $USER $PASS

# Call 3. Here you will change to the directory where you want to put or get
cd uploads/matthewklau@fas.harvard.edu_PnE9pydx/apg

# Call4.  Here you will tell FTP to put or get the file.
put $1

# End FTP Connection
bye

EOF
EOF
