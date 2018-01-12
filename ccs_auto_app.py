import pandas as pd
import xml.etree.ElementTree
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import seaborn as sns
from compute_ccs import SteppedFieldCCS

import argparse

##########################################################################
# ArgumentParser
##########################################################################
parser = argparse.ArgumentParser()

parser.add_argument(
    '--target_list_file', type=str, default='TargetList.txt',
    help='Target list file (Tab-delimited text format)')

parser.add_argument(
    '--config_file', type=str, default='config.xml',
    help='Configuration file')

parser.add_argument(
    '--data_folder', type=str, default='./',
    help='Data folder containing all the cef and meta data files')

parser.add_argument(
    '--output', type=str, default='ccs_table.tsv',
    help='Output file to save a output table')


FLAGS = {}

##########################################################################
# 
##########################################################################
def etree_to_dict(t):
    '''convert XML tree to json format
    '''
    d = {t.tag : map(etree_to_dict, t.iter())}
    d.update(('@' + k, v) for k, v in t.attrib.items())
    d['text'] = t.text
    d['tag'] = t.tag
    return d

def get_metadata(mfile, offset):
    '''read metadata file and extract the field information for each frame
    TODO: offset method (choose one frame by offset) or average in a range
    Return
        a pandas dataframe having a field information for each frame
    '''
    metadata = pd.read_csv(mfile, sep='\t')
    _list = list(metadata.drop_duplicates(subset='FrameMethodId').FrameId+offset-1)
    return metadata[metadata.FrameId.isin(_list)]

def get_target_info(target_list_file):
    '''read the target_list_file
        target_list_file: file path for a config file
    Return
        a pandas dataframe containing target information
    '''
    return pd.read_csv(target_list_file, sep='\t').fillna(method='ffill')

def get_config(config_file):
    '''read the config_file
        config_file: file path for a config file
    Return
        config_params: a dict containing a configuration information
    '''
    e = xml.etree.ElementTree.parse(config_file).getroot()
    json = etree_to_dict(e)

    config_params = {}
    adducts_params=[]
    for j in json['configuration']:
        if 'parameter' in j:
            if j['@name'] in ['mz_tolerance', 'drift_tube_length', 'old_drift_tube_length', 'neutral_mass'] :
                val = float(j['text'])
            elif j['@name'] in ['frame_offset', 'num_fields']:
                val = int(j['text'])
            elif j['@name'] != 'adducts':
                val = j['text']
            config_params[j['@name']] = val
        if 'adduct' in j:
            add = {'name':j['@name'], 'charges':j['@charges'], 'mass':float(j['@mass'])}
            adducts_params.append(add)
    config_params['adducts'] = adducts_params
    return config_params

def get_adducts(exact_mass, adducts):
    '''get the adducts mass
        exact_mass: exact mass of the target
        adducts: configuration for adducts in config_file
    Return
        adducts2mass: a dict containing information of positive and negative adducts
    '''
    adducts2mass = {'pos':{}, 'neg':{}}
    for adduct in adducts:
        charges = adduct['charges'].replace(' ','').split(',')
        for c in charges:
            charge = int(c)
            name = '[M'+c+adduct['name']+']' if abs(charge)>1 else '[M'+c[0]+adduct['name']+']'
            mass = (exact_mass + charge * adduct['mass'])/abs(charge)
            if charge > 0:
                adducts2mass['pos'][name] = mass
            elif charge < 0:
                adducts2mass['neg'][name] = mass
    return adducts2mass

def get_features_from_cef(cef, max_normalize=True):
    '''get features by reading a cef file
    '''
    e = xml.etree.ElementTree.parse(cef).getroot()
    json = etree_to_dict(e.findall('CompoundList')[0])
    idx = 0
    mppid = 0
    rst = []
    mspeaks = []
    in_ms_peaks = False
    for j in json['CompoundList']:
        if 'Compound' in j:
            mppid = j['@mppid']
        if 'Location' in j:
            mz = j['@m']
            rt = j['@rt']
            intensity = j['@y']
            dt = j['@dt']
            rst.append({'mppid':mppid, 'mz':float(mz), 'rt':float(rt), 'intensity':float(intensity), 'dt':float(dt)})
        if 'MSPeaks' in j:
            for k in j['MSPeaks']:
                if ('p' in k):
                    mspeaks.append({'mppid':mppid, 'mass':float(k['@x']), 'intensity_org':float(k['@y']), 'z':k['@z'], 's':k['@s']})
    df = pd.DataFrame(rst)
    mspeaks = pd.DataFrame(mspeaks)
    if df.shape[0] > 0:
        df['intensity_z'] = (df.intensity - df.intensity.mean())/df.intensity.std()
        if max_normalize:
            df.intensity /= df.intensity.max()
        mspeaks = mspeaks.drop_duplicates(subset='mppid', keep='first')
        df = pd.merge(mspeaks, df, left_on="mppid", right_on="mppid", how='inner')
    return df, mspeaks

def get_adducts_colors():
    return {'[M+.]':'m',
            '[M+H]':'b',
            '[M+2H]':'c',
            '[M+Na]':'r',
            '[M+K]':'g',
            '[M-H]':'y'}
        
def is_in_tolerance(x, mass, ppm):
    delta = mass * ppm * 1.0e-6
    #print(mass, delta, mass-delta, mass+delta)
    return (x >= mass - delta) & (x <= mass + delta)

def mass_error(x, mass):
    return abs(x - mass) / mass * 1e6

def find_features(features, metadata, ion_mz, ppm):
    df = features[is_in_tolerance(features.mass, ion_mz, ppm)]
    df = df.sort_values(by='intensity_z').drop_duplicates(subset='frame', keep='last')
    df = df.merge(metadata, left_on='frame', right_on='FrameMethodId', how='inner')
    df = df.sort_values(by='frame')
    return df

def get_ccs(comp_id, target_list, config_params):
    '''
    Return
        a list
    '''
    ccs_results = []

    target_info = target_list[target_list.ID==comp_id]
    rep_files = list(target_info.RawFileName)
    rep_files.sort()
    num_reps = len(rep_files)

    adducts = get_adducts(list(target_info.ExactMass)[0], config_params['adducts'])[list(target_info.Ionization)[0]]

    ##################################################
    plt.close('all')
    figs = {}
    is_filled = {}
    axis = {}
    for adduct in adducts:
        figs[adduct], axis[adduct] = plt.subplots(num_reps, sharex=True, sharey=True, figsize=(8,8))
        is_filled[adduct] = False
    figs['meta'], axis['meta'] = plt.subplots(3, sharex=True, sharey=False, figsize=(8,8))
    figs['intdist'], axis['intdist'] = plt.subplots(config_params['num_fields'], len(rep_files), sharex=True, sharey=False, figsize=(20,14))
    ##################################################

    # compute CCS for each replicate
    for r, rep_file in enumerate(rep_files):
        fname = FLAGS.data_folder + '/' + rep_file.rsplit('.', 1)[0]
        meta_file = (fname + '{0}.txt').format(config_params['suffix_meta'])
        metadata = get_metadata(meta_file, config_params['frame_offset'])
        
        ##################################################
        axis['meta'][0].plot(metadata.FrameId, metadata.ImsTemperature, label=rep_file)
        axis['meta'][0].set_ylabel('Temperature (C)')
        axis['meta'][1].plot(metadata.FrameId, metadata.ImsPressure)
        axis['meta'][1].set_ylabel('Pressure (torr)')
        axis['meta'][2].plot(metadata.FrameId, metadata.ImsField)
        axis['meta'][2].set_ylabel('E (V/cm)')
        axis['meta'][2].set_xlabel('Frame ID')
        ##################################################

        for step in range(config_params['num_fields']):
            cef_file = (fname + '{0}{1:d}.cef').format(config_params['suffix_raw'], (step+1))
            _features, _ = get_features_from_cef(cef_file)
            _features['frame'] = np.ones(_features.shape[0], dtype=np.int32) * (step+1)
            if step == 0:
                features = _features
            else:
                features = features.append(_features)
            
            ## draw m/z vs intensity
            ax = axis['intdist'][step, r]
            plot_intensity_distribution(_features, adducts, ax, config_params['mz_tolerance'])
            ax.set_xlim([0, np.log(3e6)])
        
        # compute CCS for each adducts
        for adduct in adducts:
            adduct_mass = adducts[adduct]
            ccs_features = find_features(features, metadata, adduct_mass, config_params['mz_tolerance'])
            if ccs_features.shape[0] > 0:
                ccs_property = compute(ccs_features, adduct_mass, config_params)
                print("[{0}] {1} ({2}), CCS: {3}".format(comp_id, adduct, rep_file, ccs_property['ccs']))
                tokens = comp_id.split('_')
                ccs_property['Compound_id'] = tokens[0]
                ccs_property['Ionization'] = tokens[1]
                ccs_property['adduct'] = adduct
                ccs_property['replicate'] = rep_file
                ccs_property['name'] = list(target_info.NeutralName)[0]
                ccs_property['num_features'] = ccs_features.shape[0]
                for feature in ccs_features.itertuples():
                    ccs_property['intensity_org_' + str(feature.frame)] = feature.intensity_org
                    ccs_property['intensity_z_' + str(feature.frame)] = feature.intensity_z
                    ccs_property['mass_error_' + str(feature.frame)] = mass_error(feature.mass, adduct_mass)
                ccs_results.append(ccs_property)
                ##################################################
                plot_ccs_regression_lines(axis[adduct][r], adduct, adduct_mass, ccs_features, ccs_property, title=rep_file)
                is_filled[adduct] = True
                ##################################################
    
    ##################################################                
    for adduct in adducts:
        if is_filled[adduct]:
            figs[adduct].tight_layout()
            figs[adduct].savefig(comp_id+"_"+adduct+".pdf", dpi=300)
    
    axis['meta'][0].legend()
    figs['meta'].tight_layout()
    figs['meta'].savefig(comp_id+"_meta.pdf", dpi=300)
    
    figs['intdist'].tight_layout()
    figs['intdist'].savefig(comp_id+'_intensity_dist.pdf')

    ##################################################
    return ccs_results
            
def compute(df, ion_mz, config_params):
    '''compute ccs
    '''
    params = {}
    params['temp'] = df.ImsTemperature.tolist()
    params['pressures'] = df.ImsPressure.tolist()
    # params['voltages'] = sheet_e
    # params['voltages'] = drift_tube_length * np.array(sheet_v)
    # params['arrival_time'] = sheet_dt
    # params['voltages'] = (df.ImsField*drift_tube_length).tolist()
    params['voltages'] = (df.ImsField*config_params['old_drift_tube_length']).tolist()  ## 10.869 * (78.12 / 78.236) = 10.853 for correction
    params['arrival_time'] = df.dt.tolist()
    params['neutral_mass'] = config_params['neutral_mass']
    params['drift_tube_length'] = config_params['drift_tube_length']
    params['mass'] = ion_mz
    # print(params)
    ccs, prop = SteppedFieldCCS(params=params).compute()
    # print("CCS:", ccs)
    return prop

def plot_ccs_regression_lines(axis, adduct, adduct_mass, df, prop, title, drift_tube_length=78.236):
    
    addmass = adduct_mass
    color = get_adducts_colors()[adduct]

    p_v = df.ImsPressure / (df.ImsField * drift_tube_length)
    
    p_vmax = p_v.max()
    p_vmin = p_v.min()
    axis.scatter(p_v, df.dt, c=color)
    axis.text(0.05, 0.8, '{0} {1:.6f}'.format(adduct, addmass),
            verticalalignment='bottom', horizontalalignment='left',
            transform=axis.transAxes,
            color='k', fontsize=15)
    for r in df.itertuples():
        axis.text((r.ImsPressure / (r.ImsField * drift_tube_length) + (p_vmax - p_vmin)/7), r.dt,
                  # '{0:.3f}ppm, {1:.2f}(z_score={2:.3f})'.format(mass_error(r.mass, addmass), r.intensity, r.intensity_z),
                  '{0:.3f}ppm, z_score={1:.2f}'.format(mass_error(r.mass, addmass), r.intensity_z),
                  color='k', fontsize=10)

    axis.plot(p_v, 1000 * (prop['intercept'] + prop['slope']*p_v), 'r', label='fitted line')
    axis.text(0.05, 0.65, 'r-squared:{0:.5f}'.format(prop['r_value']**2),
        verticalalignment='bottom', horizontalalignment='left',
        transform=axis.transAxes,
        color='k', fontsize=15)
    axis.text(0.05, 0.5, 'CCS:{0:.4f}'.format(prop['ccs']),
        verticalalignment='bottom', horizontalalignment='left',
        transform=axis.transAxes,
        color='k', fontsize=15)
    axis.set_title(title)
    axis.set_xlabel('Pressure/Voltages (Torr/V)')
    axis.set_ylabel('Arrival time (ms)')


def plot_intensity_distribution(features, adducts_mass, ax, ppm=50):
    colors = get_adducts_colors()
    g = sns.kdeplot(np.log(features.intensity_org), shade=True, color="b", ax=ax)
    ax.axvline(np.log(np.median(features.intensity_org)), linestyle=':')
    ax.axvline(np.log(10*np.median(features.intensity_org)), linestyle=':')
    ax.axvline(np.log(np.mean(features.intensity_org)+2*np.std(features.intensity_org)), linestyle='-.')
    for adduct in adducts_mass:
        sel = features[is_in_tolerance(features.mass, adducts_mass[adduct], ppm)]
        if sel.shape[0] > 0:
            ax.scatter(np.log(sel['intensity_org']), np.zeros(sel.shape[0]), c=colors[adduct])
    ax.set_xlabel('log(Intensity)')
    ax.set_ylabel('Density')

# def plot_mz_error_intensity(features, adducts_mass, ax, ppm=50):
#     colors = get_pos_adducts_colors()
#     ax.scatter(features.mass, features.intensity, c='k', alpha=0.2, s=3)
#     for adduct in adducts_mass:
#         sel = features[is_in_tolerance(features.mass, adducts_mass[adduct], ppm)]
#         if sel.shape[0] > 0:
#             ax.scatter(sel['mass'], sel['intensity'], c=colors[adduct])

def report(ccs_table, target_list):
    def get_stats(group):
        return {'ccs_avg': group.mean(), 'ccs_rsd': 100*group.std()/group.mean()}
    
    ccs_avg = ccs_table.groupby(['Compound_id', 'adduct'])['ccs'].apply(get_stats).unstack()
    print(ccs_avg.reset_index())
    ccs_table = pd.merge(ccs_table, ccs_avg.reset_index(), on=['Compound_id','adduct'], how='left')
    
    # save to a csv file after reordering the columns
    cols = list(ccs_table.columns)
    cols.pop(cols.index('ccs_avg'))
    cols.pop(cols.index('ccs_rsd'))
    cols.pop(cols.index('Compound_id'))
    cols.pop(cols.index('Ionization'))
    cols.pop(cols.index('adduct'))
    cols.pop(cols.index('ccs'))
    cols.pop(cols.index('mass'))
    cols.pop(cols.index('name'))

    df = ccs_table[['Compound_id','name','Ionization','adduct','mass','ccs_avg','ccs_rsd','ccs']+cols]
    df.to_csv(FLAGS.output, sep='\t')
    
def main():
    # read a list of targets
    target_list = pd.read_csv(FLAGS.target_list_file, sep='\t')
    target_list['ID']= target_list.CompID.str.cat("_"+target_list.Ionization)
    
    num_comp = np.sum(target_list.CompID.notnull())
    num_pos = np.sum(target_list.Ionization=='pos')
    num_neg = np.sum(target_list.Ionization=='neg')

    target_list = target_list.fillna(method='ffill')

    # read a set of configuration parameters
    config_params = get_config(FLAGS.config_file)

    # compounds
    compound_ids = list(target_list.ID.drop_duplicates())
    assert len(compound_ids) == num_pos+num_neg,\
        "Please check if there are duplicates in CompID and its Ionization"

    cids = list(target_list.CompID.drop_duplicates())
    print('Number of compounds: {0} (pos:{1},neg:{2})'.format(len(cids), num_pos, num_neg))
    
    ccs_results = []
    for comp_id in compound_ids:
        # compute ccs for this compound
        ccs_results += get_ccs(comp_id, target_list, config_params)
    ccs_table = pd.DataFrame(ccs_results)
    report(ccs_table, target_list)

if __name__ == '__main__':
    FLAGS = parser.parse_args()
    main()
