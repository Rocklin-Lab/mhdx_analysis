import pandas as pd
import numpy as np
import os
import argparse
import json
import jax.numpy as jnp


def load_data(f_hx, f_hb, param_table_path, cooperativity_dict_path):
    df_hx = pd.read_json(f_hx)
    df_hbonds = pd.read_json(f_hb)
    param_table = pd.read_json(param_table_path)

    with open(cooperativity_dict_path, 'r') as json_file:
        cooperativity_dict = json.load(json_file)

    return df_hx, df_hbonds, param_table, cooperativity_dict


def preprocess_data(df_hx, df_hbonds):
    df = pd.merge(df_hx, df_hbonds[["sequence", "n_hb_bb_all", "PF"]], on="sequence", how="left")
    df["free_energy"] = df["free_energy"].apply(lambda x: sorted(x, reverse=True))
    df["n_exch_res"] = df["free_energy"].apply(len)
    df["free_energy_integrated"] = df["free_energy"].apply(np.sum)
    df["free_energy_integrated_per_res"] = df["free_energy_integrated"] / df["n_exch_res"]
    df["free_energy_integrated_measurable_per_hb"] = df["free_energy_integrated"] / df["n_hb_bb_all"]
    df["ratio_n_hb_n_exch"] = df["n_hb_bb_all"] / df["n_exch_res"]
    df["ratio_n_measurable_obs_rates_n_exch"] = df["n_measurable_obs_rates"] / df["n_exch_res"]
    df['ratio_n_measurable_intrinsic_rates_n_exch'] = df['n_measurable_intrinsic_rates'] / df['n_exch_res']
    return df


def apply_param_table_to_data(df, param_table, query):
    model = param_table.model.iloc[0]
    fit_model_terms = {}

    for pf in param_table['PF'].unique():
        params = param_table.query('PF == @pf').iloc[0]
        fit_model_terms[pf] = {
            'offset': params['offset'],
            'scaling_exp': params['scaling_exp'],
            'scaling_factor_dg': params['scaling_factor_dg'],
            'scaling_factor_nc': params['scaling_factor_nc'],
            'hb_exp': params['hb_exp']
        }

    for pf in df['PF'].unique():
        if pf not in fit_model_terms:
            fit_model_terms[pf] = fit_model_terms['all']

    def calculate_fit_pred_pf(row):
        pf_terms = fit_model_terms[row['PF']]
        return float(
            pf_terms['scaling_factor_dg'] * jnp.power(row['dg_mean'] - pf_terms['offset'], pf_terms['scaling_exp']) *
            jnp.power(row['ratio_n_hb_n_exch'], pf_terms['hb_exp']) + pf_terms['scaling_factor_nc'] * row['netcharge']
        )

    def calculate_fit_pred(row):
        all_terms = fit_model_terms['all']
        return float(
            all_terms['scaling_factor_dg'] * jnp.power(row['dg_mean'] - all_terms['offset'], all_terms['scaling_exp']) *
            jnp.power(row['ratio_n_hb_n_exch'], all_terms['hb_exp']) + all_terms['scaling_factor_nc'] * row['netcharge']
        )

    df[f'{model}_pred_pf'] = df.apply(calculate_fit_pred_pf, axis=1)
    df[f'{model}_pred'] = df.apply(calculate_fit_pred, axis=1)
    df[f'cooperativity_model_pf'] = df['free_energy_integrated_per_res'] - df[f'{model}_pred_pf']
    df[f'cooperativity_model_global'] = df['free_energy_integrated_per_res'] - df[f'{model}_pred']

    return df.query(query).reset_index(drop=True)


def normalize_cooperativity(df, cooperativity_dict):
    pf_column = 'cooperativity_model_pf'
    normalized_pf_column = 'normalized_cooperativity_model_pf'

    df[normalized_pf_column] = df.apply(
        lambda x: (
            (x[pf_column] - cooperativity_dict[x['PF']]['mean']) / cooperativity_dict[x['PF']]['std']
            if x['PF'] in cooperativity_dict and cooperativity_dict[x['PF']]['std'] != 0
            else x[pf_column]
        ),
        axis=1
    )

    global_column = 'cooperativity_model_global'
    normalized_global_column = 'normalized_cooperativity_model_global'
    global_mean = cooperativity_dict['global']['mean']
    global_std = cooperativity_dict['global']['std']

    df[normalized_global_column] = df[global_column].apply(
        lambda x: (x - global_mean) / global_std if global_std != 0 else x
    )

    return df


def save_output(df, output_path):
    if not os.path.isdir(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    df.to_json(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Normalized Cooperativities")
    parser.add_argument('--hx', required=True, help="Path to deduplicated.json")
    parser.add_argument('--hb', required=True, help="Path to hbonds.json")
    parser.add_argument('--param_table', required=True, help="Path to param_table.json")
    parser.add_argument('--cooperativity_dict', required=True, help="Path to cooperativity_std_mean_dict.json")
    parser.add_argument('--output', default="results/df_cooperativity.json", help="Path to save the final dataframe")
    parser.add_argument('--query',
                        default="group in ['group_1: measurable unmerged','group_2: measurable merged'] & dg_mean > 2 & ((-3 < normalized_cooperativity_model_global < 3) | (-3 < normalized_cooperativity_model_pf < 3))",
                        help="Query string to filter the dataframe")

    args = parser.parse_args()

    df_hx, df_hbonds, param_table, cooperativity_dict = load_data(args.hx, args.hb, args.param_table,
                                                                  args.cooperativity_dict)
    df = preprocess_data(df_hx, df_hbonds)
    df = apply_param_table_to_data(df, param_table, args.query)
    df = normalize_cooperativity(df, cooperativity_dict)
    save_output(df.query(args.query), args.output)

    print(f"Final dataframe saved to {args.output}")
