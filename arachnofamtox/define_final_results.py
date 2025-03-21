import operator
from operator import itemgetter


def pred(dic, maxormin):
    # Returns min or max value of a dictionary
    res = "results"
    maxormin = str(maxormin)
    if maxormin == "max":
        res = max(dic.items(), key=operator.itemgetter(1))[0]
    elif maxormin == "min":
        res = min(dic.items(), key=operator.itemgetter(1))[0]
    return res


def define_final(famname):
    # Define final name of family
    final = famname
    final_name = {
        "KTx": "Scorpion_KTx",
        "theraphotoxin": "theraphotoxin",
        "NaTx": "Scorpion_NaTx",
        "Metalloprotease": "Metalloprotease",
        "Metalloproteinase": "Metalloprotease",
        "ctenitoxin": "ctenitoxin",
        "hexatoxin": "hexatoxin",
        "oxotoxin": "oxotoxin",
        "zodatoxin": "zodatoxin",
        "agatoxin": "agatoxin",
        "segestritoxin": "segestritoxin",
        "filistatoxin": "filistatoxin",
        "lycotoxin": "lycotoxin",
        "sparatoxin": "sparatoxin",
    }

    for i in final_name.keys():
        if i in famname:
            final = final_name[i]
    return final


def gen_final_class_and_evalue(dic, method):
    Target_toxin = {
        "beta-": "Voltage-Gated Sodium Channels (Nav)",
        "gamma-": "Voltage-Gated Potassium Channels (Kv)",
        "delta-": "Voltage-Gated Sodium Channels (Nav)",
        "kappa-": "Voltage-Gated Potassium Channels (Kv)",
        "mu-": "Voltage-Gated Sodium Channels (Nav)",
        "tau-": "Transient Receptor Potential Vanilloid 1 (TRPV1)",
        "omega-": "Voltage-Gated Calcium Channels (Cav)",
        "M-": "Membranolytic Activity",
        "alpha-latrotoxin": "Presynaptic Neurons",
        "Holocyclotoxin": "Presynaptic Neurons",
    }

    Target = {
        "KTx": "Voltage-Gated Potassium Channels (Kv)",
        "NaTx": "Voltage-Gated Sodium Channels (Nav)",
        "CAP": "Antigen-5, Pathogen-related, including ion channels",
        "Lectin": "Glycoproteins on cell surfaces",
        "PeptidaseS1": "Serine Proteases",
        "PhospholipaseA2": "Cell membranes",
        "Cystatin": "Inhibits cysteine proteases",
        "Kunitz-type": "Protease Inhibitors, Neurotoxic Effects",
        "TickDefensin": "Microbial cell membranes",
        "NDBP": "Various ion channels and receptors",
        "AstacinLikeMetalloprotease": "Extracellular matrix components",
        "TickMetalloprotease": "Extracellular Matrix",
        "Cytolytic": "Cell Membranes",
        "Dermonecrotic": "Sphingomyelin Phosphodiesterase D",
        "GlycosylHydrolase_Hyaluronidase": "Hyaluronic acid in the extracellular matrix",
        "Lipocalin": "Small hydrophobic molecules",
        "Prokineticin": "Unknown",
        "Holocyclotoxin": "Presynaptic Neurons",
        "TickMetalloprotease": "Extracellular Matrix",
        "VenomMetalloproteinase": "Extracellular Matrix",
    }

    target = "-"

    if method == "hmm":  # for hmm get min value
        result_class = pred(dic, "min")
    elif method == "pssm":  # for pssm get max value
        result_class = pred(dic, "max")

    possible_class = str(result_class).lower()
    if "toxin" in possible_class:
        for i in Target_toxin.keys():
            possible_target = i.lower()
            if possible_target in possible_class:
                target = Target_toxin[i]  # define target for neurotoxins
    else:
        for i in Target.keys():
            possible_target = i.lower()
            if possible_target in possible_class:
                target = Target[i]  # define target for other toxins

    result_class_final = define_final(result_class)  # final name of family
    result_evalue_final = dic[result_class]  # final e-value
    return result_class, result_class_final, result_evalue_final, target


def gen_score(result_hmm_final, result_pssm_final):
    # Generate score from PSSM and HMM results
    result_final = "None"
    score = 0
    if result_hmm_final == result_pssm_final:
        result_final = result_hmm_final
        score = 5
    elif result_hmm_final == "None":
        result_final = result_pssm_final
        score = 3
    elif result_pssm_final == "None":
        result_final = result_hmm_final
        score = 2
    else:
        result_final = result_pssm_final  # use PSSM as final
        score = 4
    return result_final, score


def filter_results(result_final, score):
    write = True
    if result_final == "GlycosylHydrolase_Hyaluronidase" and score == 2:
        write = False
    elif result_final == "alpha-latrotoxin" and score == 2:
        write = False
    elif result_final == "Dermonecrotic" and score == 2:
        write = False
    return write


def define_general_class(result_final):
    VP = [
        "Metalloprotease",
        "CAP",
        "Cystatin",
        "GlycosylHydrolase_Hyaluronidase",
        "Kunitz-type",
        "Lectin",
        "NDBP",
        "PeptidaseS1",
        "PhospholipaseA2",
        "Prokineticin",
        "TickDefensin",
        "TickLipocalin",
    ]

    TX = [
        "alpha-latrotoxin",
        "Scorpion_KTx",
        "Scorpion_NaTx",
        "theraphotoxin",
        "ctenitoxin",
        "hexatoxin",
        "oxotoxin",
        "zodatoxin",
        "agatoxin",
        "segestritoxin",
        "filistatoxin",
        "lycotoxin",
        "sparatoxin",
    ]

    if result_final in VP:
        gen_class = "Venom Protein"
    else:
        gen_class = "Toxin"
    return gen_class


def define_target_final(target_hmm, target_pssm):
    target_final = "-"
    if target_hmm == target_pssm:
        target_final = target_pssm
    elif target_pssm == "-" and target_hmm != "-":
        target_final = target_hmm
    elif target_hmm == "-" and target_pssm != "-":
        target_final = target_pssm
    else:
        target_final = target_pssm
    return target_final
