import pandas as pd
import xml.etree.ElementTree
import os
import numpy as np
import matplotlib.pyplot as plt
import glob
import seaborn as sns
from compute_ccs import SteppedFieldCCS
from shortest_for_ccs import get_possible_ccs_values
import time

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

parser.add_argument(
    '--r2_threshold', type=float, default=0.99,
    help='threshold value for r2')

parser.add_argument(
    '--num_isotopes_threshold', type=int, default=2,
    help='threshold value for num_isotopes')

parser.add_argument(
    '--intensity_rank_threshold', type=int, default=3,
    help='threshold value for peak intensity rank in m/z window')

parser.add_argument(
    '--maxint', action='store_true',
    help='select max intensive peaks for ccs computation')

parser.add_argument(
    '--format', type=str, default='cef',
    help='file format for the features, e.g., cef or mzmine')

FLAGS = {}

##########################################################################
# 
##########################################################################

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

def get_metadata(mfile, offset, ax=None, label=None):
    '''read metadata file and extract the field information for each frame
    TODO: offset method (choose one frame by offset) or average in a range
    Return
        a pandas dataframe having a field information for each frame
    '''
    try:
        metadata = pd.read_csv(mfile, sep='\t')
        _list = list(metadata.drop_duplicates(subset='FrameMethodId').FrameId+offset-1)
        filtered = metadata[metadata.FrameId.isin(_list)]
        ##################################################
        if ax is not None:
            ax[0].plot(metadata.FrameId, metadata.ImsTemperature, label=label)
            ax[0].scatter(filtered.FrameId, filtered.ImsTemperature, label=None)
            ax[0].set_ylabel('Temperature (C)')
            ax[1].plot(metadata.FrameId, metadata.ImsPressure)
            ax[1].scatter(filtered.FrameId, filtered.ImsPressure)
            ax[1].set_ylabel('Pressure (torr)')
            ax[2].plot(metadata.FrameId, metadata.ImsField)
            ax[2].scatter(filtered.FrameId, filtered.ImsField)
            ax[2].set_ylabel('E (V/cm)')
            ax[2].set_xlabel('Frame ID')
        ##################################################
        return filtered
    except Exception as e:
        return None


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
        # num_isotopes = mspeaks.mppid.value_counts()
        mspeaks['num_isotopes'] = mspeaks.groupby('mppid')['mppid'].transform('count')
        # print(num_isotopes)
        mspeaks = mspeaks.drop_duplicates(subset='mppid', keep='first')
        # mspeaks['num_isotopes'] = num_isotopes
        
        df = pd.merge(mspeaks, df, left_on="mppid", right_on="mppid", how='inner')
    return df, mspeaks#, num_isotopes


def get_features_from_mzmine_csv(csv, max_normalize=True):
    df = pd.read_csv(csv)
    if df.shape[0] == 0: return df, None
    col_id = [c for c in df.columns if c.endswith('row ID')][0]
    col_area = [c for c in df.columns if c.endswith('Peak area')][0]
    col_height = [c for c in df.columns if c.endswith('Peak height')][0]
    col_mz = [c for c in df.columns if c.endswith('Peak m/z')][0]
    col_dt = [c for c in df.columns if c.endswith('Peak RT')][0]
    if 'calibrated_ccs' in df.columns:
        cols = [col_id, col_mz, col_dt, col_area, col_height, 'calibrated_ccs']
        colnames = ['mppid', 'mass', 'dt', 'intensity', 'height', 'calibrated_ccs']
    else:
        cols = [col_id, col_mz, col_dt, col_area, col_height]
        colnames = ['mppid', 'mass', 'dt', 'intensity', 'height']
    
    df = df[cols].copy()
    df.columns = colnames
    df['intensity_org'] = df.intensity
    if df.shape[0] > 0:
        df['intensity_z'] = (df.intensity - df.intensity.mean())/df.intensity.std()
        if max_normalize:
            df.intensity /= df.intensity.max()
    return df, None


def get_features(file, max_normalize=True, format='cef'):
    if format=='cef': return get_features_from_cef(file, max_normalize)
    elif format=='mzmine': return get_features_from_mzmine_csv(file, max_normalize)
    else: print('File format: {0}. This tool doesn\'t support this file format.')
    return None, None

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

def find_features_maxint(features, metadata, ion_mz, ppm):
    df = features[is_in_tolerance(features.mass, ion_mz, ppm)]
    if df.shape[0] == 0: return df
    
    #  if 'frame' column in metadata, delete it
    if 'frame' in metadata.columns: del metadata['frame']

    df = df.sort_values(by='intensity_z').drop_duplicates(subset='frame', keep='last')
    df = df.merge(metadata, left_on='frame', right_on='FrameMethodId', how='inner')
    df = df.sort_values(by='frame')
    return df

def find_features(features, metadata, ion_mz, ppm,
                  threshold_num_isotopes=2,
                  threshold_intensity_rank=3):
    if 'num_isotopes' in features.columns:
        df = features[is_in_tolerance(features.mass, ion_mz, ppm) & (features.num_isotopes>=threshold_num_isotopes)]
    else:
        df = features[is_in_tolerance(features.mass, ion_mz, ppm)]
    if df.shape[0] == 0: return df
    
    # filter out small peaks by ranking threshold
    rankings = df.groupby('frame')['intensity_org'].rank(ascending=False)
    df = df[rankings<=threshold_intensity_rank]

    # for f in frames_too_many_features:
    #     filter_by_intensity_rank(df, f, threshold_intensity_rank)

    #  if 'frame' column in metadata, delete it
    if 'frame' in metadata.columns: del metadata['frame']

    # df = df.sort_values(by='intensity_z').drop_duplicates(subset='frame', keep='last')
    df = df.merge(metadata, left_on='frame', right_on='FrameMethodId', how='inner')
    # df = df.sort_values(by='frame')
    # df.to_csv("test_{0:.5f}.txt".format(ion_mz),sep="\t")
    return df

def filter_by_intensity_rank(df, frame, threshold_intensity_rank=3):
    temp = df[df.frame == frame]
    # print(df)
    # print(frame, temp.intensity_org)
    np.argsort(temp.intensity_org)

def ccs_filter(ccs_list):
    # remove the redundant regression lines which share the same start nodes(features)
    first_peaks = []
    last_peaks = []
    for ccs in ccs_list:
        first_peaks.append(int(ccs.mppid[0]))
        last_peaks.append(int(ccs.mppid[-1]))
    
    ufirst_peaks = list(np.unique(first_peaks))
    ulast_peaks = list(np.unique(last_peaks))
    if len(ufirst_peaks) < len(ccs_list):
        print("len(ufirst_peaks) < len(ccs_list)", len(ufirst_peaks),len(ccs_list))
        _ccs_list = []
        for u in ufirst_peaks:
            idx_list = np.where(first_peaks == u)[0]
            if idx_list.shape[0] > 1:
                best_r2 = 0
                best_ccs_u = None
                for ii in idx_list:
                    if (best_r2 < ccs_list[ii].r2):
                        best_ccs_u = ccs_list[ii]
                        best_r2 = ccs_list[ii].r2
                if best_ccs_u != None:
                    _ccs_list.append(best_ccs_u)
            else:
                _ccs_list.append(ccs_list[idx_list[0]])
        return _ccs_list

    elif len(ulast_peaks) < len(ccs_list):
        print("len(ulast_peaks) < len(ccs_list)", len(ulast_peaks),len(ccs_list))
        print("ulast_peaks", ulast_peaks)
        print("last_peaks", last_peaks)
        _ccs_list = []
        for u in ulast_peaks:
            idx_list = np.where(last_peaks == u)[0]
            print('idx_list',u, idx_list)
            if idx_list.shape[0] > 1:
                best_r2 = 0
                best_ccs_u = None
                for ii in idx_list:
                    if (best_r2 < ccs_list[ii].r2):
                        best_ccs_u = ccs_list[ii]
                        best_r2 = ccs_list[ii].r2
                if best_ccs_u != None:
                    _ccs_list.append(best_ccs_u)
            else:
                _ccs_list.append(ccs_list[idx_list[0]])
        return _ccs_list
    else:
        return ccs_list
        
    # find the ccs values of earlist molecules
    pass


def files_not_enough(fname, config_params, format='cef'):
    meta_file = (fname + '{0}.txt').format(config_params['suffix_meta'])
    if not os.path.isfile(meta_file):
        print("[ERROR] a metadata file doesn't exist:", meta_file)
        return True
    for step in range(config_params['num_fields']):
        if format=='cef': ffile = (fname + '{0}{1:d}.cef').format(config_params['suffix_raw'], (step+1))
        else: ffile = (fname + '{0}{1:d}.mzML_c_dc_de.csv').format(config_params['suffix_raw'], (step+1))
        if not os.path.isfile(ffile):
            print("[ERROR] a feature file doesn't exist:", ffile)
            return True
    return False


def get_ccs(comp_id, target_list, config_params,
            threshold_r2,
            threshold_num_isotopes,
            threshold_intensity_rank):
    '''
    Return
        a list
    '''
    ccs_results = []

    time_for_feature_finding = 0

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
    figs['intdist'], axis['intdist'] = plt.subplots(config_params['num_fields'], num_reps, sharex=True, sharey=False, figsize=(20,14))
    ##################################################

    # compute CCS for each replicate
    try:
        for r, rep_file in enumerate(rep_files):
            fname = FLAGS.data_folder + '/' + rep_file.rsplit('.', 1)[0]
            
            if files_not_enough(fname, config_params, FLAGS.format):
                ccs_prop = dict()
                tokens = comp_id.rsplit('_', 1)
                ccs_prop['Compound_id'] = tokens[0]
                ccs_prop['Ionization'] = tokens[1]
                # ccs_prop['adduct'] = adduct
                ccs_prop['replicate'] = rep_file
                ccs_prop['name'] = list(target_info.NeutralName)[0]
                ccs_prop['CAS'] = list(target_info.CAS)[0]
                ccs_prop['comments'] = "couldn't find some files to compute CCS"
                ccs_results.append(ccs_prop)
                continue
            
            meta_file = (fname + '{0}.txt').format(config_params['suffix_meta'])
            metadata = get_metadata(meta_file, config_params['frame_offset'], ax=axis['meta'], label=rep_file)
            
            # collecting features
            features = []
            for step in range(config_params['num_fields']):
                if FLAGS.format=='cef': ffile = (fname + '{0}{1:d}.cef').format(config_params['suffix_raw'], (step+1))
                else: ffile = (fname + '{0}{1:d}.mzML_c_dc_de.csv').format(config_params['suffix_raw'], (step+1))
                
                _features, _ = get_features(ffile, format=FLAGS.format)
                if _features.shape[0] > 0:
                    _features['frame'] = np.ones(_features.shape[0], dtype=np.int32) * (step+1)
                    features.append(_features)
                    
                    ## draw m/z vs intensity
                    if num_reps == 1:
                        ax = axis['intdist'][step]
                    else:
                        ax = axis['intdist'][step, r]
                    plot_intensity_distribution(_features, adducts, ax, config_params['mz_tolerance'])
                    ax.set_xlim([0, np.log(3e6)])
                else:
                    print("[ERROR] This file has no features: {0}".format(ffile))
            
            if len(features) == 0: continue
            features = pd.concat(features)

            # compute CCS for each adducts
            # print(features)
            # print("features size:", features.shape)

            for adduct in adducts:
                adduct_mass = adducts[adduct]
                start_time = time.time()
                
                if (FLAGS.maxint):
                    ccs_features_within_mz = find_features_maxint(features, metadata, adduct_mass, config_params['mz_tolerance'])
                else:
                    ccs_features_within_mz = find_features(features, metadata, adduct_mass, config_params['mz_tolerance'],
                                                       threshold_num_isotopes=threshold_num_isotopes,
                                                       threshold_intensity_rank=threshold_intensity_rank)
                if ccs_features_within_mz.shape[0] > 0:
                    ccs_list = get_possible_ccs_values(ccs_features_within_mz,
                                                       adduct_mass,
                                                       old_drift_tube_length=config_params['old_drift_tube_length'],
                                                       drift_tube_length=config_params['drift_tube_length'],
                                                       neutral_mass=config_params['neutral_mass'],
                                                       threshold_n_fields=3,
                                                       threshold_r2=threshold_r2)
                    # filtering should be done based on ccs values of across all 3 replicates
                    # Note: i am not sure if r2 is a good metric to do this.
                    ccs_list = ccs_filter(ccs_list)

                    if len(ccs_list) > 0:
                        tokens = comp_id.rsplit('_', 1)
                        for ccs in ccs_list:
                            ccs_prop = ccs.to_dict()
                            print("[{0}] {1} ({2}), CCS: {3}({4})".format(comp_id, adduct, rep_file, ccs_prop['ccs'], ccs_prop['r2']))
                            ccs_prop['Compound_id'] = tokens[0]
                            ccs_prop['Ionization'] = tokens[1]
                            ccs_prop['adduct'] = adduct
                            ccs_prop['replicate'] = rep_file
                            ccs_prop['name'] = list(target_info.NeutralName)[0]
                            ccs_prop['CAS'] = list(target_info.CAS)[0]
                            ccs_results.append(ccs_prop)

                        if num_reps == 1:
                            _tmp_ax = axis[adduct]
                        else:
                            _tmp_ax = axis[adduct][r]

                        ##################################################
                        plot_ccs_regression_lines2(
                            _tmp_ax,
                            adduct,
                            adduct_mass,
                            ccs_features_within_mz,
                            ccs_list,
                            title=rep_file,
                            drift_tube_length=config_params['drift_tube_length'])
                        is_filled[adduct] = True
                        ##################################################
                # ccs_features = find_features_maxint(features, metadata, adduct_mass, config_params['mz_tolerance'])
                # time_for_feature_finding += (time.time() - start_time)
                # if ccs_features.shape[0] > 0:
                #     ccs_property = compute(ccs_features, adduct_mass, config_params)
                #     print("[{0}] {1} ({2}), CCS: {3}".format(comp_id, adduct, rep_file, ccs_property['ccs']))
                #     tokens = comp_id.rsplit('_', 1)
                #     ccs_property['Compound_id'] = tokens[0]
                #     ccs_property['Ionization'] = tokens[1]
                #     ccs_property['adduct'] = adduct
                #     ccs_property['replicate'] = rep_file
                #     ccs_property['name'] = list(target_info.NeutralName)[0]
                #     ccs_property['CAS'] = list(target_info.CAS)[0]
                #     ccs_property['num_features'] = ccs_features.shape[0]
                #     for feature in ccs_features.itertuples():
                #         ccs_property['intensity_org_' + str(feature.frame)] = feature.intensity_org
                #         ccs_property['intensity_z_' + str(feature.frame)] = feature.intensity_z
                #         ccs_property['mass_error_' + str(feature.frame)] = mass_error(feature.mass, adduct_mass)
                #     ccs_results.append(ccs_property)
                #     ##################################################
                #     if num_reps == 1:
                #         plot_ccs_regression_lines(axis[adduct], adduct, adduct_mass, ccs_features, ccs_property, title=rep_file, drift_tube_length=config_params['drift_tube_length'])
                #     else:
                #         plot_ccs_regression_lines(axis[adduct][r], adduct, adduct_mass, ccs_features, ccs_property, title=rep_file, drift_tube_length=config_params['drift_tube_length'])
                #     is_filled[adduct] = True
                #     ##################################################
        
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
    except Exception as e:
        print ("ERROR: {0} ({1})".format(e.strerror, rep_file))
    # print('Total time for feature finding: {0} sec/compound(e.g., 3 reps and 6 adducts)'.format(time_for_feature_finding))
    return ccs_results
            
def compute(df, ion_mz, config_params):
    '''compute ccs
    '''
    params = {}
    params['temp'] = df.ImsTemperature.tolist()
    params['pressures'] = df.ImsPressure.tolist()
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

# def plot_ccs_regression_lines2(axis, adduct, adduct_mass, df, prop, title, drift_tube_length=78.236):
def plot_ccs_regression_lines2(
                    axis,
                    adduct,
                    adduct_mass,
                    df,
                    ccs_list,
                    title,
                    drift_tube_length):
    addmass = adduct_mass
    color = get_adducts_colors()[adduct]

    p_v = df.ImsPressure / (df.ImsField * drift_tube_length)
    
    p_vmax = p_v.max()
    p_vmin = p_v.min()

    for r in df.itertuples():
        axis.scatter(r.ImsPressure / (r.ImsField * drift_tube_length), r.dt,
            # c=color, s=np.log(r.intensity_org), alpha=0.2)
            c=color, s=1000*r.intensity, alpha=0.2)

    axis.text(0.05, 0.8, '{0} {1:.6f}'.format(adduct, addmass),
            verticalalignment='bottom', horizontalalignment='left',
            transform=axis.transAxes,
            color='k', fontsize=15)

    for ccs in ccs_list:
        prop = ccs.to_dict()
        pv = [ccs.pressures[i] / (ccs.fields[i] * drift_tube_length) for i in range(len(ccs.pressures))]
        dt_diff = [abs(ccs.arrival_time[i-1]-ccs.arrival_time[i]) for i in range(1,len(ccs.arrival_time))]
        for i, f in enumerate(ccs.fields):
            axis.text((pv[i] + (p_vmax - p_vmin)/7), ccs.arrival_time[i],
                      '{0:.3f}ppm, z_score={1:.2f}'.format(ccs.mass_ppm_error[i], ccs.intensity_z[i]),
                      color='k', fontsize=10)
            # axis.scatter(pv[i], ccs.arrival_time[i], s=np.log(ccs.intensity_org[i]), c=color)
            axis.scatter(pv[i], ccs.arrival_time[i], s=1000*ccs.intensity[i], c=color, alpha=0.8)
        
        axis.text(min(pv)-2*(p_vmax - p_vmin)/7, min(ccs.arrival_time)-0.25*min(dt_diff),
            'CCS:{0:.4f}(r2:{1:.5f})'.format(prop['ccs'], prop['r2']),
            color='r', fontsize=10)

        axis.plot(p_v, 1000 * (prop['intercept'] + prop['slope']*p_v), 'r', label='fitted line')
    axis.set_title(title)
    axis.set_xlabel('Pressure/Voltages (Torr/V)')
    axis.set_ylabel('Arrival time (ms)')


def plot_intensity_distribution(features, adducts_mass, ax, ppm=50):
    if features.shape[0] > 0:
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


def report(ccs_table, target_list):
    if ccs_table.shape[0] == 0:
        print("Unfortunately, we couldn't find any good CCS values.")
        return
    def get_stats(group):
        return {'ccs_avg': group.mean(), 'ccs_rsd': 100*group.std()/group.mean(), 'ccs_count': group.count()}
    
    ccs_avg = ccs_table.groupby(['Compound_id', 'adduct'])['ccs'].apply(get_stats).unstack()
    print(ccs_avg.reset_index())
    ccs_table = pd.merge(ccs_table, ccs_avg.reset_index(), on=['Compound_id','adduct'], how='left')
    
    print(ccs_table.head())
    # save to a csv file after reordering the columns
    cols = list(ccs_table.columns)
    if 'ccs_avg' in cols:
        cols.pop(cols.index('ccs_avg'))
    else:
        ccs_table['ccs_avg'] = np.nan
    if 'ccs_rsd' in cols:
        cols.pop(cols.index('ccs_rsd'))
    else:
        ccs_table['ccs_rsd'] = np.nan
    cols.pop(cols.index('Compound_id'))
    cols.pop(cols.index('Ionization'))
    cols.pop(cols.index('adduct'))
    cols.pop(cols.index('ccs'))
    cols.pop(cols.index('adduct_mass'))
    cols.pop(cols.index('name'))
    newcols = ['Compound_id','name','Ionization','adduct','adduct_mass','ccs_avg','ccs_rsd','ccs']+cols
    
    df = ccs_table[newcols]
    # df = ccs_table
    df.to_csv(FLAGS.output, sep='\t')
    
def main():
    # read a list of targets
    target_list = pd.read_csv(FLAGS.target_list_file, sep='\t')
    target_list['ID']= target_list.CompID.str.cat("_"+target_list.Ionization)
    
    target_list = target_list.fillna(method='ffill')

    ## e.g., S00001.b if you have a same compound id but different versions.
    # num_comp = list(pd.DataFrame(target_list.CompID.str.split('\.').tolist(), columns = ['CompID','ver']).CompID.drop_duplicates())
    num_comp = list(target_list.CompID.drop_duplicates())
    num_pos = np.sum(pd.DataFrame(target_list.ID.drop_duplicates().str.rsplit('_',1).tolist(), columns = ['id','ion']).ion=='pos')
    num_neg = np.sum(pd.DataFrame(target_list.ID.drop_duplicates().str.rsplit('_',1).tolist(), columns = ['id','ion']).ion=='neg')

    # read a set of configuration parameters
    config_params = get_config(FLAGS.config_file)

    # compounds
    compound_ids = list(target_list.ID.drop_duplicates())
    assert len(compound_ids) == num_pos+num_neg,\
        "Please check if there are duplicates in CompID and its Ionization"

    print('Number of compounds: {0} (pos:{1},neg:{2})'.format(len(num_comp), num_pos, num_neg))
    
    ccs_results = []
    start_time = time.time()
    for comp_id in compound_ids:
        # compute ccs for this compound
        ccs_results += get_ccs(comp_id, target_list, config_params,
                               threshold_r2=FLAGS.r2_threshold,
                               threshold_num_isotopes=FLAGS.num_isotopes_threshold,
                               threshold_intensity_rank=FLAGS.intensity_rank_threshold)
        print('[{0}] {1} sec'.format(comp_id, (time.time()-start_time)))
    print('Total time: {0} sec/compound(e.g., 3 reps)'.format((time.time()-start_time)/len(compound_ids)))
    ccs_table = pd.DataFrame(ccs_results)
    report(ccs_table, target_list)

if __name__ == '__main__':
    FLAGS = parser.parse_args()
    print("options:", FLAGS)
    main()
