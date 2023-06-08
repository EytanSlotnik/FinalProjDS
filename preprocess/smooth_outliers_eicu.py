import pandas as pd


def load_bp(num_of_pats=300):
    drug = pd.read_csv("../data/eICU/infusiondrug.csv")
    diagnosis = pd.read_csv("../data/eICU/diagnosis.csv")
    patients_weight = pd.read_csv("../data/eICU/patient.csv")[["patientunitstayid", "admissionweight"]]
    patient = pd.read_csv("../data/eICU/patient.csv")
    drug = drug.merge(patients_weight, on="patientunitstayid")

    # add hospitalID and unitType and wardID to the drug table corelated to patientunitstayid
    drug["hospitalid"] = drug["patientunitstayid"].map(patient.set_index("patientunitstayid")["hospitalid"])
    drug["unittype"] = drug["patientunitstayid"].map(patient.set_index("patientunitstayid")["unittype"])
    drug["wardid"] = drug["patientunitstayid"].map(patient.set_index("patientunitstayid")["wardid"])

    # filter drug from NA in drugname, patientunitstayid, hospitalid, unittype, wardid
    drug = drug[~drug['drugname'].isna()]
    drug = drug[~drug['patientunitstayid'].isna()]
    drug = drug[~drug['hospitalid'].isna()]
    drug = drug[~drug['unittype'].isna()]
    drug = drug[~drug['wardid'].isna()]

    # filter drugs with other drug names
    other_drugs = drug[drug['drugname'].str.startswith(("Epinephrine", "Dopamine", "Vesopressin"))]

    # filter drug to include only drugnames that start with "Norepinephrine"
    drug = drug[drug['drugname'].str.startswith("Norepinephrine")]

    # remove patients who recieved other vasopressors
    drug = drug[~drug['patientunitstayid'].isin(other_drugs['patientunitstayid'])]

    # sort drug by patientunitstayid and drugstartoffset
    drug = drug.sort_values(by=["patientunitstayid", "infusionoffset"])

    bp = pd.read_csv("../preprocess/filtered_bp_eicu.csv")
    bp = bp.sort_values(by=["stay_id", "cur_bp_time"])
    bp = bp[bp["stay_id"].isin(drug["patientunitstayid"])]
    bp["next_bp_time"] = bp.groupby("stay_id")["cur_bp_time"].shift(-1)
    bp["interval"] = bp["next_bp_time"] - bp["cur_bp_time"]
    drug["drugrate"] = pd.to_numeric(drug["drugrate"], errors='coerce')
    bp["age"] = pd.to_numeric(bp["stay_id"].map(patient.set_index("patientunitstayid")["age"]), errors='coerce')
    bp["hospitalid"] = bp["stay_id"].map(patient.set_index("patientunitstayid")["hospitalid"])
    bp["unittype"] = bp["stay_id"].map(patient.set_index("patientunitstayid")["unittype"])
    bp["wardid"] = bp["stay_id"].map(patient.set_index("patientunitstayid")["wardid"])

    # find the number of patients per hospital
    bp_hosp = bp.groupby(["hospitalid"]).agg({"stay_id": ["nunique"]}).sort_values(by=("stay_id", "nunique"))

    bp['num_of_pats'] = bp['hospitalid'].map(bp_hosp[('stay_id', 'nunique')])

    # bp_big will consist of patients from hospitalid with more then 400 differnt patientsunitstayid
    bp_big = bp[bp["hospitalid"].isin(drug.groupby("hospitalid")["patientunitstayid"].nunique().sort_values()[
                                          drug.groupby("hospitalid")[
                                              "patientunitstayid"].nunique().sort_values() > num_of_pats].index)]
    return bp_big


def smooth_outliers(big_bp: pd.DataFrame, threshold_constant=3):
    for i in range(3, 4):
        big_bp["rolling_" + str(i)] = big_bp.groupby("stay_id")["cur_bp"].rolling(i).mean().reset_index(0, drop=True)
        # replace NaN with the original value
        big_bp["rolling_" + str(i)] = big_bp["rolling_" + str(i)].fillna(big_bp["cur_bp"])
        big_bp["rolling_" + str(i) + "_res"] = big_bp["cur_bp"] - big_bp["rolling_" + str(i)]
        big_bp["rolling_" + str(i) + "_res"] = big_bp["rolling_" + str(i) + "_res"].abs()

        # add column with the median of the residuals
        big_bp["rolling_" + str(i) + "_res_median"] = big_bp.groupby("stay_id")["rolling_" + str(i) + "_res"].transform(
            "median")
        # add column with outliers removed (outliers are defined as residuals that are 10 times the median)
        big_bp["rolling_" + str(i) + "_res_median"] = big_bp["rolling_" + str(i) + "_res_median"] * threshold_constant
        big_bp["smooth_" + str(i)] = big_bp["cur_bp"]
        big_bp.loc[big_bp["rolling_" + str(i) + "_res"] > big_bp["rolling_" + str(i) + "_res_median"], "smooth_" + str(
            i)] = big_bp["rolling_" + str(i)]
    return big_bp


if __name__ == "__main__":
    big_bp = load_bp(num_of_pats=50)
    big_bp.to_csv("../preprocess/big_bp_eicu.csv", index=False)
    smooth_bp = smooth_outliers(big_bp)
    smooth_bp.to_csv("../preprocess/smooth_bp_eicu.csv", index=False)
