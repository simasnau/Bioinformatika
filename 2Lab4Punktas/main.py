from ete3 import Tree

camel_virus = 'MT085175.1'
generated = Tree('treeOutput.txt', format=1)

print("Prieš:")
print(generated)


generated.set_outgroup(camel_virus)

print("Po:")
print(generated)