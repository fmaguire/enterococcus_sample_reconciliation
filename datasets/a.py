import lxml.etree as etree

x = etree.parse('PRJNA604849_sra_data.xml')
print(etree.tostring(x, pretty_print=True))
