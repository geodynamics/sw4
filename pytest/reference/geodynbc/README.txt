# Approach to making test for whi file and geodynbc
Requires a new test script verify_geodynbc.py
---
1. Create station files at whi nodes using loh1-h100-mr-whi.in (from scec and hdf5). Station files will go in loh1-h100-mr-whi/
2. Convert station files to whi file using whiwriter.py. File is named loh1-h100-mr.whi
3. Use resulting whi file as input on geodynbc input line in loh-h100-mr-geodynbc.in and create station files. Station files will go in loh-h100-mr-geodynbc/
4. Check station files against reference (a la hdf5). Will create verify_geodynbc.py
