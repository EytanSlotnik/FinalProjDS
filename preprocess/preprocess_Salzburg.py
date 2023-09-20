import pandas as pd
import gzip
import struct
import os
import csv



current_dir = os.getcwd()
os.chdir(current_dir)


def SetRawValues(row, dictwriter, n):
    t = int(row['Offset'])
    data = bytes.fromhex(row["rawdata"][2:])  # deserialize hex string to bytes
    for i in range(int(len(data) / 4)):
        if data[i * 4] == 0 and data[i * 4 + 1] == 0 and data[i * 4 + 2] == 0 and data[i * 4 + 3] == 0:
            continue  # no null values
        n = n + 1  # new primary key
        newrow = row.copy()
        if newrow["DataID"] != '703':
            continue
        del newrow["DataID"]
        del newrow["rawdata"]  # not needed
        del newrow["cnt"]  # not needed
        del newrow["id"]  # not needed
        # newrow["id"] = n  # primary key
        newrow["Val"] = struct.unpack('<f', data[i * 4:i * 4 + 4])[0]  # bytes to float
        newrow["Offset"] = t + i * 60  # new offset
        dictwriter.writerow(newrow)
    return n


n = 0
with open('data_float_m_703.csv', 'w', newline='') as csvfile:
    dict_writer = csv.DictWriter(csvfile, ['id', 'CaseID', 'DataID', 'Offset', 'Val'])
    dict_writer = csv.DictWriter(csvfile, ['CaseID', 'Offset', 'Val'])
    dict_writer.writeheader()
# counter = 10000
    with gzip.open('data_float_h.csv.gz', 'rt') as gzf:
        for row in csv.DictReader(gzf):
#             # counter -= 1
#             # if counter == 0:
#             #     break
            n = SetRawValues(row, dict_writer, n)
            if (n % 1000 == 0):
                 print("Processing entry " + str(n))

d703 = pd.read_csv('data_float_m_703.csv')
d703 = d703[d703['Offset'] % 300 == 0]
d703 = d703[d703['Offset'] <= 172800]

d_cases = pd.read_csv('cases.csv.gz', compression='gzip')
d_cases['WeightOnAdmission'] /= 1000

medication = pd.read_csv('medication.csv.gz', compression='gzip')


d1562 = medication[medication['DrugID'] == 1562]
d1562 = d1562[d1562['OffsetDrugEnd'] <= 172800]

id_weight_dict = dict(zip(d_cases['CaseID'], d_cases['WeightOnAdmission']))

intersection_ids = list(set(d703['CaseID'].unique()).intersection(d1562['CaseID'].unique()))


# Create a boolean mask to filter rows from d703 and d1562
mask1 = d703['CaseID'].isin(intersection_ids)
mask2 = d1562['CaseID'].isin(intersection_ids)

# Apply the boolean mask to filter rows from d703 and d1562
df1 = d703.loc[mask1, ['CaseID', 'Offset', 'Val']]
df2 = d1562.loc[mask2, ['CaseID', 'Offset', 'OffsetDrugEnd', 'AmountPerMinute']]

# Concatenate df1 and df2
df3 = pd.concat([df1, df2], axis=0)

# Concatenate df3 with df4
df3.reset_index(drop=True, inplace=True)

df3['Offset'] /= 60
df3['OffsetDrugEnd'] /= 60
df3['AmountPerMinute'] *= 1000000

df3.rename(columns={
    'CaseID' : 'stay_id',
    'Val': 'cur_bp',       # Rename 'Val' to 'Value'
    'Offset': 'cur_bp_time',  # Rename 'Drug Rate' to 'Rate'
    'AmountPerMinute': 'drugrate'  # Rename 'Offset' to 'Offset Range'
}, inplace=True)


death_dict = {2202: 'ICU survivor', 2212: 'Unknown', 2215: 'ICU non-survivor'}
death_per_pat_dict = dict(zip(d_cases['CaseID'], d_cases['DischargeState']))


df3['DischargeState'] = df3['stay_id'].map(death_per_pat_dict)
df3['DischargeState'] = df3['DischargeState'].map(death_dict)


df3['weight'] = df3['stay_id'].map(id_weight_dict)


df3.to_csv('five_minutes_resolution.csv', index=False)

