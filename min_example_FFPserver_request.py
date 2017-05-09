#!/usr/bin/env python
# Copyright (C) 2015, 2016 OpenEye Scientific Software
import os,sys,json,requests

baseurl = "http://130.180.63.34:8089"

print( "Get all databases loaded" )
response = requests.get( baseurl )
data = response.json()

for entry in data[ 'databases' ]:
    print( entry )

#print( "\nGet more information about the loaded databases" )
#for entry in data['databases']:
#    url = "%s/%s" %( baseurl, entry )
#    response = requests.get( url )
#    res = response.json()
#    print( res )

print( "\nPerform a query on first db" )
url = "%s/%s/hitlist?smiles=%s&oformat=csv" %( baseurl, data['databases'][0],"c1ccccc1CCC1CCCCC1")
response = requests.get( url )
print(response.content)
#for entry in res:
#    print( "MolID: %-20s  SimVal: %.4f" %( entry, res[ entry ] ) )

#print( "\nPerform a query on first db with cutoff and maxhits" )
#url = "%s/%s/hitlist?smiles=%s&cutoff=%f&maxhits=%d" %( baseurl, data['databases'][0],"c1ccccc1CCC1CCCCC1", 0.3, 5)
#response = requests.get( url )
#res = response.json()
#for entry in res:
#    print( "MolID: %-20s  SimVal: %.4f" %( entry, res[ entry ] ) )
