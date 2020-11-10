from operon_analyzer import rules, visualize
import sys

COLOR_SELECTOR = { "cas1": "lightblue",
                    "cas2": "seagreen",
                    "cas3": "gold",
                    "cas4": "springgreen",
                    "cas5": "darkred",
                    "cas6": "thistle",
                    "cas7": "coral",
                    "cas8": "red",
                    "cas9": "palegreen",
                    "cas10": "yellow",
                    "cas11": "tan",
                    "cas12": "orange",
                    "cas13": "saddlebrown",
                    "tnsA": "navy",
                    "tnsB": "blue",
                    "tnsC": "purple",
                    "tnsD": "white",
                    "tniQ": "teal",
                    "tnsE": "gray",
                    "glmS": "orchid",
                    "CRISPR array": "blue",
                    "casphi": "olive"
                    }

operon_list = []
operons = visualize.build_operon_dictionary(sys.stdin)
for operon in operons.values():
    operon_list.append(operon)
visualize.plot_operons(operons=operon_list, output_directory=sys.argv[1], feature_colors=COLOR_SELECTOR, nucl_per_line=25000)