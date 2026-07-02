import os, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from matplotlib.lines import Line2D

BASE='./data/Supplementary_Data_S2.xlsx'
OUT='./outputs/AUC_regenerated_figures'
os.makedirs(OUT, exist_ok=True)
allp=pd.read_excel(BASE, sheet_name='02_All_Pairs')
top=pd.read_excel(BASE, sheet_name='03_Top_Ranked')
casewise=pd.read_excel(BASE, sheet_name='06_Casewise')
curves=pd.read_excel(BASE, sheet_name='07_Dose_Curves')

# Ensure numeric
for df in [allp, top, casewise, curves]:
    for col in df.columns:
        if col not in ['Case','Drug_ID','Drug_Name_Short','Drug_Name_Full','Directional_Status','Direction','TopHit_Directional_Discordance','Large_Rank_Shift_ge50']:
            try: df[col]=pd.to_numeric(df[col])
            except Exception: pass

# S9
rho_all, p_all=spearmanr(allp['Delta_Min'], allp['Delta_AUC'])
rho_top, p_top=spearmanr(top['Delta_Min'], top['Delta_AUC'])
fig, axes=plt.subplots(2,2, figsize=(13.5,10.5), dpi=200)
ax=axes[0,0]
ax.scatter(allp['Delta_Min'], allp['Delta_AUC'], s=12, alpha=0.45, edgecolor='none')
ax.axhline(0,color='0.78',lw=1)
ax.axvline(0,color='0.78',lw=1)
ax.plot([-105,105],[-105,105],'--',color='0.35',lw=1.2)
ax.set_xlim(-105,105); ax.set_ylim(-105,105)
ax.set_xlabel('ΔMin (GMI − PBS)')
ax.set_ylabel('ΔAUC (GMI − PBS)')
ax.set_title(f'All compound–case pairs\nSpearman ρ={rho_all:.3f}, P={p_all:.1e}', fontsize=12)
ax.text(-0.15,1.08,'A',transform=ax.transAxes,fontsize=18,fontweight='bold', va='top')

ax=axes[0,1]
colors=np.where(top['TopHit_Directional_Discordance'].astype(bool), '#d62728', '#1f77b4')
ax.scatter(top['Delta_Min'], top['Delta_AUC'], s=24, alpha=0.8, c=colors, edgecolor='white', linewidth=0.25)
ax.axhline(0,color='0.78',lw=1)
ax.axvline(0,color='0.78',lw=1)
ax.plot([-105,105],[-105,105],'--',color='0.35',lw=1.2)
ax.set_xlim(-105,105); ax.set_ylim(-105,105)
ax.set_xlabel('ΔMin (top-ranked)')
ax.set_ylabel('ΔAUC')
ax.set_title(f'Top-ranked ΔMin pairs\nSpearman ρ={rho_top:.3f}, P={p_top:.1e}', fontsize=12)
legend_elements=[
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#1f77b4', markersize=7, label='ΔAUC concordant with ΔMin'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#d62728', markersize=7, label='ΔAUC discordant with ΔMin')]
ax.legend(handles=legend_elements, frameon=False, fontsize=9, loc='lower right')
ax.text(-0.15,1.08,'B',transform=ax.transAxes,fontsize=18,fontweight='bold', va='top')

ax=axes[1,0]
rankdiff=allp['Rank_Diff_Sensitized'].dropna()
ax.hist(rankdiff, bins=45, color='#1f77b4', alpha=0.85)
ax.axvline(0,color='0.35',lw=1.2)
ax.axvline(-50,color='0.5',lw=1.0,ls='--')
ax.axvline(50,color='0.5',lw=1.0,ls='--')
med=np.nanmedian(rankdiff); q1=np.nanpercentile(rankdiff,25); q3=np.nanpercentile(rankdiff,75)
ax.text(0.98,0.93, f'Median = {med:.0f}\nIQR = {q1:.0f} to {q3:.0f}', transform=ax.transAxes, ha='right', va='top', fontsize=10)
ax.set_xlabel('Rank difference (AUC rank − ΔMin rank)')
ax.set_ylabel('Compound–case count')
ax.set_title('Distribution of rank differences', fontsize=12)
ax.text(-0.15,1.08,'C',transform=ax.transAxes,fontsize=18,fontweight='bold', va='top')

ax=axes[1,1]
casewise=casewise.sort_values('Case')
bars=ax.bar(casewise['Case'], casewise['Spearman_rho_all'], color='#1f77b4')
ax.set_ylim(0,1.0)
ax.set_ylabel('Spearman ρ')
ax.set_xlabel('Case ID')
ax.set_title('Case-wise ΔMin–ΔAUC agreement', fontsize=12)
ax.tick_params(axis='x', rotation=45)
for b,v in zip(bars,casewise['Spearman_rho_all']):
    ax.text(b.get_x()+b.get_width()/2, v+0.02, f'{v:.2f}', ha='center', va='bottom', fontsize=8)
ax.text(-0.15,1.08,'D',transform=ax.transAxes,fontsize=18,fontweight='bold', va='top')

for ax in axes.flat:
    ax.grid(False)
    ax.spines[['top','right']].set_visible(True)
fig.tight_layout(rect=[0,0,1,0.98], h_pad=2.2, w_pad=2.2)
fig.savefig(os.path.join(OUT,'Supplementary_Figure_S9.png'), dpi=600, bbox_inches='tight')
fig.savefig(os.path.join(OUT,'Supplementary_Figure_S9.pdf'), bbox_inches='tight')
plt.close(fig)

# S10: chosen representative examples (explicit drug IDs to resolve duplicates)
examples=[
    ('A127','D176','Moboc'),
    ('A127','D173','Gilte'),
    ('A051','D140','Cedir'),
    ('A166','D116','Creno'),
    ('A174','D116','Creno'),
    ('A174','D067','Etopo'),
    ('A166','D035','Dacar'),
    ('A127','D160','Eribu'),
    ('A127','D167','Tepot'),
    ('A173','D006','Entre'),
    ('A043','D118','Praci'),
    ('A173','D008','Masit'),
]
fig, axes=plt.subplots(3,4, figsize=(15,10), dpi=200, sharex=True, sharey=True)
for i,(case,drugid,short) in enumerate(examples):
    ax=axes.flat[i]
    m=curves[(curves['Case']==case)&(curves['Drug_ID']==drugid)].sort_values('Dose_Position_Low_to_High')
    row=allp[(allp['Case']==case)&(allp['Drug_ID']==drugid)].iloc[0]
    x=m['Dose_Position_Low_to_High'].to_numpy()
    ax.plot(x, m['PBS'].to_numpy(), marker='o', lw=1.6, ms=4.5, label='PBS')
    ax.plot(x, m['GMI'].to_numpy(), marker='o', lw=1.6, ms=4.5, label='GMI')
    ax.set_ylim(-5,105); ax.set_xlim(0.7,8.3)
    ax.set_xticks(range(1,9))
    ax.set_title(f'{i+1}. {case} {short}\nΔMin={row.Delta_Min:.1f}, ΔAUC={row.Delta_AUC:.1f}', fontsize=9)
    if i % 4 == 0:
        ax.set_ylabel('Normalized viability (%)')
    if i >= 8:
        ax.set_xlabel('Relative dose level\n(low → high)')
    ax.tick_params(axis='both', labelsize=8)

handles, labels=axes.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.995), ncol=2, frameon=False, fontsize=11)
fig.suptitle('Representative discordant examples: normalized viability across relative dose levels', y=0.965, fontsize=14, fontweight='bold')
fig.text(0.01, 0.012, 'ΔMin and ΔAUC are defined as (GMI − PBS). Negative values indicate increased sensitivity with GMI; positive values indicate attenuation with GMI.', fontsize=9)
fig.tight_layout(rect=[0.02,0.04,1,0.935], h_pad=2.0, w_pad=1.4)
fig.savefig(os.path.join(OUT,'Supplementary_Figure_S10.png'), dpi=600, bbox_inches='tight')
fig.savefig(os.path.join(OUT,'Supplementary_Figure_S10.pdf'), bbox_inches='tight')
plt.close(fig)

print('Saved files:')
for f in sorted(os.listdir(OUT)):
    print(os.path.join(OUT,f))
