from ete3 import Tree

camel_virus = 'lcl|Query_3460_4901-8458_MN514967.1_Dromedary_camel_coronavirus_HKU23_isolate_DcCoV-HKU23/camel/Nigeria/NV1385/2016'
generated = Tree('treeOutput.txt', format=1)

print("Prieš:")
print(generated)


generated.set_outgroup(camel_virus)

print("Po:")
print(generated)

generated.write(format=1, outfile="modifiedTree.txt")
generated.render("modifiedTree.png")