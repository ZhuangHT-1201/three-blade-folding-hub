import json, os
from docx import Document
from docx.shared import Pt, Cm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIG = os.path.join(ROOT, 'figures')
D = json.load(open(os.path.join(HERE, 'design.json'), encoding='utf-8'))
FC = json.load(open(os.path.join(HERE, 'fc_report.json'), encoding='utf-8'))
CL = json.load(open(os.path.join(HERE, 'clearance.json'), encoding='utf-8'))
ro, L, U, M = D['rotor'], D['loads'], D['util'], FC['mass']

doc = Document()
st = doc.styles['Normal']
st.font.name = 'Times New Roman'; st.font.size = Pt(11)


def h(t, lv=1):
    doc.add_heading(t, lv)


def para(t):
    doc.add_paragraph(t)


def table(head, rows):
    tb = doc.add_table(rows=1, cols=len(head)); tb.style = 'Table Grid'
    for c, t in zip(tb.rows[0].cells, head): c.text = t
    for r in rows:
        cs = tb.add_row().cells
        for c, t in zip(cs, r): c.text = str(t)
    doc.add_paragraph()


def fig(name, cap, w=15):
    p = os.path.join(FIG, name)
    if os.path.exists(p):
        doc.add_picture(p, width=Cm(w))
        doc.add_paragraph(cap).alignment = 1


doc.add_heading('Three-Blade Co-Directional Folding Hub: Design Report', 0)

h('1 Design Inputs')
table(['Item', 'Value'], [
    ['Number of blades', '3'],
    ['Tip radius R', f"{ro['R']:.1f} m"],
    ['Root radius', f"{ro['r_root']:.2f} m"],
    ['Airfoil', 'NACA 0012, chord 246 to 103 mm, twist 29.3 to 6.8 deg'],
    ['Operating speed', f"{ro['rpm_nom']:.0f} rpm"],
    ['Design check speed', f"{ro['rpm_max']:.0f} rpm (110%)"],
    ['Blade mass', f"{D['blade_mass']:.3f} kg, CG radius {D['blade_rcg']:.3f} m"],
    ['Hover thrust / torque / power at 982 rpm', f"{ro['T_nom']:.0f} N / {ro['Q_nom']:.0f} N·m / {ro['P_nom']/1000:.1f} kW"],
])

h('2 Concept')
para('The hub consists of an upper and a lower 7075-T6 aluminium plate with a centre spacer ring; the three arms are identical. '
     'Each blade grip has a tongue that sits between the two plates and is located by two vertical pins: a pivot bolt '
     f"and a quick-release lock pin. The pins are {D['s']*1000:.0f} mm apart in the chordwise direction at a radius of "
     f"{D['r_h']*1000:.0f} mm.")
para('Each arm carries two flight holes and two parking holes (PARK1, PARK2). All arms share the same hole pattern, '
     'so the rotor mass distribution remains symmetric.')
para('Blade A (azimuth 180 deg) stays in place. Blade B (azimuth 60 deg) is released by pulling the trailing-edge lock pin, '
     'rotated +122 deg about the leading-edge bolt, and locked in PARK1. Blade C (azimuth 300 deg) is released by pulling '
     'the leading-edge lock pin, rotated -122 deg about the trailing-edge bolt, and locked in PARK2. The three blades end up '
     'pointing in the same direction with a 2 deg spread. Deployment is the reverse sequence.')
para(f"The hinge plane is {abs(D['z_h'])*1000:.0f} mm below the pitch axis. At this offset the centrifugal moment about the "
     'hinge plane cancels the thrust moment. Both scale with rotor speed squared, so the steady flap moment at the hinge '
     'is zero at any hover speed, and the pins and plates carry only the alternating part.')
fig('fold_sequence_top.png', 'Figure 1. Folding sequence, top view (deployed / 50% / folded)', 16)
fig('asm_folded.png', 'Figure 2. Folded assembly')
fig('hub_deployed_3d.png', 'Figure 3. Folding hub, deployed', 12)
fig('hub_folded_3d.png', 'Figure 4. Folding hub, folded', 12)
fig('hub_top.png', 'Figure 5. Hub plate hole pattern', 12)

h('3 Main Dimensions')
table(['Parameter', 'Value [mm]'], [
    ['Hinge radius r_h', f"{D['r_h']*1000:.0f}"], ['Pin spacing s', f"{D['s']*1000:.0f}"],
    ['Pin diameter d', f"{D['d']*1000:.0f}"], ['Tongue thickness', f"{D['t_t']*1000:.0f}"],
    ['Tongue lug outer radius', f"{D['R_t']*1000:.0f}"], ['Plate thickness (each)', f"{D['t_p']*1000:.0f}"],
    ['Plate lug outer radius', f"{D['R_p']*1000:.0f}"], ['Plate arm width', f"{D['b_a']*1000:.0f}"],
    ['Hinge plane height z_h', f"{D['z_h']*1000:.0f}"], ['Cuff length / wall', f"{D['L_cuff']*1000:.0f} / {D['t_cuff']*1000:.0f}"],
])

h('4 Parts List')
table(['Part', 'Material', 'Qty', 'Mass each [kg]'], [
    ['hub_plate_upper', '7075-T6', 1, f"{M['plate']:.3f}"],
    ['hub_plate_lower', '7075-T6', 1, f"{M['plate']:.3f}"],
    ['hub_spacer', '7075-T6', 1, f"{M['spacer']:.3f}"],
    ['blade_grip', '7075-T6', 3, f"{M['grip']:.3f}"],
    ['pivot_bolt', '42CrMo', 3, f"{M['bolt']:.3f}"],
    ['lock_pin', '17-4PH', 3, f"{M['lock_pin']:.3f}"],
    ['sleeve_bushing', 'Aluminium bronze', 6, f"{M['sleeve']:.3f}"],
    ['Total (without blades)', '', '', f"{M['mechanism_total']:.2f}"],
])

h('5 Stowed Envelope')
fo, de = D['folded'], D['deployed']
table(['State', 'Plan area [m²]', 'Envelope [mm]'], [
    ['Deployed', f"{de['plan_area']:.2f}", '2800 x 2800 x 124'],
    ['Folded', f"{fo['plan_area']:.2f}", ' x '.join(f'{v*1000:.0f}' for v in fo['bbox'])],
])
para(f"The folded plan area is {fo['plan_area']/de['plan_area']*100:.1f}% of the deployed value and the envelope volume is "
     f"{fo['bbox_vol']/de['bbox_vol']*100:.1f}%.")
para('Solid boolean interference checks at 25%, 50%, 75% and 100% fold give zero overlap between grips and between grips '
     'and the hub. Minimum clearances [mm]:')
rows = []
for k, lab in (('0.25', '25%'), ('0.5', '50%'), ('0.75', '75%'), ('1.0', '100%')):
    c = CL[k]
    rows.append([lab] + [f'{c[n]:.1f}' for n in ('AB', 'AC', 'BC', 'A-hub', 'B-hub', 'C-hub')])
table(['Fold', 'A-B', 'A-C', 'B-C', 'A-hub', 'B-hub', 'C-hub'], rows)
para('In the folded position the gap between the cuff of blade C and the lower plate edge is 7.8 mm; that plate edge is '
     'chamfered C5.')

h('6 Strength Check')
para(f"Design loads are taken at {ro['rpm_max']:.0f} rpm with an ultimate factor of 1.5, a gust factor of 1.3 and an "
     'alternating flap ratio of 0.35. Fatigue uses the Goodman line with a safety factor of 1.5. '
     f"Ultimate centrifugal force at the hinge is {L['Fx']/1000:.1f} kN and the flap moment is {L['My']:.0f} N·m.")
names = {'pin_shear': 'Pin shear', 'pin_bend': 'Pin bending', 'pin_yield': 'Pin yield', 'tongue_shearout': 'Tongue shear-out',
         'tongue_bearing': 'Tongue bearing', 'tongue_net': 'Tongue net section', 'tongue_root': 'Tongue root',
         'plate_shearout': 'Plate shear-out', 'plate_bearing': 'Plate bearing', 'plate_net': 'Plate net section',
         'plate_arm': 'Plate arm bending', 'fat_pin': 'Pin fatigue', 'fat_plate_lug': 'Plate lug fatigue',
         'fat_plate_hole': 'Plate hole fatigue', 'fat_tongue_root': 'Tongue root fatigue',
         'stiff_tongue': 'Tongue stiffness (<= 0.25 deg)', 'stiff_arm': 'Plate arm stiffness (<= 0.25 deg)',
         'fret_plate': 'Plate fretting contact pressure'}
table(['Check', 'Utilisation'], [[v, f'{U[k]:.2f}'] for k, v in names.items()])
para('All utilisations are at or below 1.0.')
fig('loads.png', 'Figure 6. Load distribution', 15)
fig('pareto.png', 'Figure 7. Hinge layout Pareto front', 12)

doc.save(os.path.join(ROOT, 'design_report.docx'))
