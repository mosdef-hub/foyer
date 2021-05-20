"""Utility to modify a silica/oplsaa xml."""
from lxml import etree

root = etree.fromstring(open("oplsaa_with_silica.xml", "r").read())
atomtypes = root[0]
bondforce = root[1]
anglforce = root[2]
rbtorsion = root[3]
periodictorsion = root[4]
nonbonded = root[5]
all_atomtypes = []
for xml_line in atomtypes:
    all_atomtypes.append(xml_line.attrib["name"])

all_atomtypes = set(all_atomtypes)

# Modify overrides to only include valid overrides
for xml_line in atomtypes:
    overrides = xml_line.attrib.get("overrides", "")
    overrides_list = overrides.split(",")
    remove_overrides = []
    new_overrides = []
    if len(overrides) > 0:
        new_overrides = [a for a in overrides_list if a in all_atomtypes]
    xml_line.attrib["overrides"] = ",".join(a for a in new_overrides)


# Trim rbtorsions
to_remove = []
for xml_line in rbtorsion:
    possible_types = [a for a in xml_line.attrib if "type" in a]
    torsion_types = {val for k, val in xml_line.attrib.items() if "type" in k}
    # torsion_types = {xml_line.attrib['type1'], xml_line.attrib['type2'],
    # xml_line.attrib['type3'], xml_line.attrib['type4']}
    if not torsion_types.issubset(all_atomtypes):
        to_remove.append(xml_line)
for remove_this in to_remove:
    rbtorsion.remove(remove_this)
out_string = etree.tostring(root)
with open("output.xml", "wb") as f:
    f.write(out_string)
