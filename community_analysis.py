def com_cover_patient(com, patient_genes):
    if len(set(com).intersection(set(patient_genes))) == len(com):
        return True
    return False
