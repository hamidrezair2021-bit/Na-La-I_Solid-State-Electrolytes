import ternary
from matplotlib import pyplot as plt
from pymatgen.core import Composition
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

elements = ['I','Na','La']

with open('VASPRUN_SOURCE_hull_report.txt','r') as hull:
    convex_hull_data=hull.read().split('\n')

Convex_hull=[]
get_data=False
for i in convex_hull_data:
    if 'Entry' in i:
        get_data = True
    if get_data:
        Convex_hull.append(i.split())

Convex_hull_energy0 = sorted([[entry[0],float(entry[1])] for entry in Convex_hull[1:-1]],key=lambda x:x[1])
Convex_hull_energy1=[]
for compound in Convex_hull_energy0:
    formula=compound[0].split('_')[0]
    comp=Composition(formula).get_el_amt_dict()
    comp_keys=comp.keys()
    comp1=[]
    for atom in elements:
        if atom in comp_keys:
            comp1.append(comp[atom])
        else:
            comp1.append(0)
    Convex_hull_energy1.append([formula,comp1,compound[1]])
Convex_hull_energy=[]
compound = []
for conv in Convex_hull_energy1:
    if conv[0] in compound:
        pass
    else:
        norm=sum(conv[1])
        Convex_hull_energy.append([conv[0],[100*i/norm for i in conv[1]],conv[2]])
        compound.append(conv[0])

ref=[]
for compound in Convex_hull_energy:
    if compound[2]==0.0:
        if 100.0 in compound[1]:
            pass
        else:
            ref.append(compound)
print(ref)

scale = 100
fig, tax = ternary.figure(scale=scale)
fig.set_size_inches(16, 11)  # Sets width and height in inches
#fig.set_dpi(600)
tax.boundary(linewidth=2)
tax.gridlines(color="blue", multiple=10)

#tax.set_title("Convex Hull", fontsize=14)
#tax.left_axis_label()
#tax.right_axis_label()
#tax.bottom_axis_label()
#tax.ticks(axis='lbr', multiple=100, linewidth=1)
tax.get_axes().set_frame_on(False)
tax.clear_matplotlib_ticks()
tax.set_axis_limits({'b': [0, 100], 'l': [0, 100], 'r': [0, 100]})
tax.get_ticks_from_axis_limits()



#tax.set_custom_ticks(fontsize=18, offset=0.01, linewidth=0)


#plt.scatter(x, y, c=energies, cmap='RdYlGn_r', vmin=0, vmax=150)

def formulas(compound):
    compo=Composition(compound).get_el_amt_dict()
    compo_keys=compo.keys()
    
    stio=r''
    for i in compo_keys:
        stio+=i
        if compo[i]==1:
            pass
        else:
            stio+='$_'+str(round(compo[i]))+'$'
    return stio
custom_ticks = {'l': {}, 'b': {}, 'r': {}}
for i in ref:
    # Extract your identifier and coordinates from the reference tuple
    label_name = formulas(i[0])
    coord_a, coord_b, coord_c = i[1][0], i[1][1], i[1][2]

    if coord_a == 0:
        # Map the position coordinate to its custom text string
        custom_ticks['l'][coord_c] = label_name
    elif coord_b == 0:
        print(coord_a)
        #tax._ticks['b'] = [elements[0] if x == coord_a else "" for x in tax._ticks['b']]

    else:
        custom_ticks['r'][coord_b] = label_name
tax._ticks['b'] = [elements[0] if x == 100 else "" for x in tax._ticks['b']]
tax._ticks['l'] = [elements[2] if x == 100 else "" for x in tax._ticks['l']]
tax._ticks['r'] = [elements[1] if x == 100 else "" for x in tax._ticks['r']]
tax.set_custom_ticks(fontsize=18, offset=0.01, linewidth=0)

for i in Convex_hull_energy:
    tax.scatter([i[1]],c=i[2],marker='o' if i[2]==0 else 'D',cmap='RdYlGn_r', vmin=0, vmax=Convex_hull_energy[-1][2]*0.35)

for i in range(len(ref)):
    for j in range(i+1,len(ref)):
        if j!=len(ref):
            if ref[i][1].index(0.0)!=ref[j][1].index(0.0):
                point1=ref[i][1]
                point2=ref[j][1]
                tax.plot([point1,point2],label=formulas(ref[i][0])+'+'+formulas(ref[j][0]))
tax.legend(prop={'size': 18},loc='upper left')
norm = Normalize(vmin=0, vmax=1000*Convex_hull_energy[-1][2]*0.35)
sm = ScalarMappable(norm=norm, cmap='RdYlGn_r')
sm.set_array([]) # Keeping the data buffer blank prevents plotting artifacts

# 4. Draw the colorbar using the virtual mappable 'sm' linked to the axis
fig = plt.gcf()
ax = tax.get_axes() # Extracts the core matplotlib geometry from the ternary container
cbar = fig.colorbar(sm, ax=ax, pad=0.08, shrink=0.8)
cbar.set_label('Energy above Hull (meV)', rotation=270, labelpad=20, fontsize = 18)
cbar.ax.tick_params(labelsize=18)

print(Convex_hull_energy)

plt.savefig('Convex_hull_energy.png')
plt.show()



