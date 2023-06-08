import pandas as pd

if __name__ == '__main__':
    bp = pd.read_csv("../preprocess/filtered_bp_eicu.csv")
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

    # save drug table
    drug.to_csv("../data/eICU/filtered_drug.csv", index=False)

    # read filtered drug table
    # drug = pd.read_csv("../data/eICU/filtered_drug.csv")
