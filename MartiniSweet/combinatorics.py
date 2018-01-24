import numpy as np
import itertools
#~ from mechanize import Browser
#~ import urllib
#~ import urllib2
#~ import webbrowser
#~ import requests
#~ import mechanize


encoding={}
with open('Hexose') as f:
   for line in f:
      key,value=line.strip().split(":")
      encoding[key]=value
      
x=encoding.values()
combinations=[list(p) for p in itertools.product(x, repeat=3)]

# 1. We will use pyranose rings only
# 2. we only use 1-4 linkage
# alpha + beta
# L and D
isomers=["L","D"] 
ring_form=["p","f"] #["pyranose","furanose"]
linkages=["1-4","1-3","1-2","1-6"]
anomers=["a","b"] # alpha and beta
sugar1="Man"
#~ conf1="D"
#~ ring="p"
n=0
tri_combinations=[list(p) for p in itertools.product(x, repeat=3)]
tri_isomers=[list(i) for i in itertools.product(isomers, repeat=3)]
tri_ringtypes=[list(i) for i in itertools.product(ring_form, repeat=3)]
di_linkages=[list(i) for i in itertools.product(linkages, repeat=2)]
tri_anomers=[list(i) for i in itertools.product(anomers, repeat=3)]
for c in itertools.product(tri_isomers,tri_combinations,tri_ringtypes,tri_anomers,di_linkages):
   # now we zip it together to get correct names
   monomers=[''.join(collection) for collection in zip(*c[:-1])]
   links=c[-1]
   name=monomers[0]+links[0]+monomers[1]+links[1]+monomers[2]+'1-OH'
   #~ print name
   
   n=n+1
   #~ if n>9:
      #~ break
print n
quit()
# now drawing: first we decide which sugars we connect:
#~ for c in combinations:
   # Then we decide for each component what is the isomer
   #~ print c
   #~ for sugar in c:
      #~ print 
   #~ for link in linkages:
      #~ for link in linkages:
         #~ for conf in configurations:
            #~ i+=1
#~ print i
#~ request=sugar1+conf1+ring+"1"+"-OH"
#~ DManpb1-OH
#~ url = "http://glycam.org/url?glycam="+request
#~ br = mechanize.Browser()
#~ br.set_handle_robots(False) # ignore robots
#~ br.open(url)
#~ br.follow_link("http://glycam.org/tools/molecular-dynamics/oligosaccharide-builder/download-files")
#~ br.select_form(name="theForm")


#~ data = urllib.urlencode({'q': 'Python'})
#~ results = urllib2.urlopen(url, data)
#~ with open("results.html", "w") as f:
    #~ f.write(results.read())

#~ webbrowser.open("results.html")
#~ print [i for i in results if 'Submit' in i]