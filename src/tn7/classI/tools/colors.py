# Okabe Ito color palette
black = (0.0, 0.0, 0.0)
orange = (0.9019607843137255, 0.6235294117647059, 0.0)
lightblue = (0.33725490196078434, 0.7058823529411765, 0.9137254901960784)
green = (0.0, 0.6196078431372549, 0.45098039215686275)
yellow = (0.9411764705882353, 0.8941176470588236, 0.25882352941176473)
blue = (0.0, 0.4470588235294118, 0.6980392156862745)
red = (0.8352941176470589, 0.3686274509803922, 0.0)
purple = (0.8, 0.4745098039215686, 0.6549019607843137)
# And an extra gray
gray = (0.2, 0.2, 0.2)


feature_colors = {'transposase|transposition|transposon|integrase|integration|resolvase|recombinase|recombination|IS\d+|(T|t)np': blue,
                  'cas(1$|2|4)': green,
                  'cas(3|9|10|12|13)': red,
                  'cas6|case|csy4': red,
                  'cas5|casd|csc1|csy2|csf3|csm4|csx10|cmr3': lightblue,
                  'cas7|casc|csd2|csc2|csy3|csf2|csm3|csm5|cmr1|cmr6|cmr4': purple,
                  'cas8|casa|csh1|csd1|cse1|csy1|csf1': yellow,
                  'cas11|casb|cse2|csm2|cmr5': red,
                  '(CRISPR )?array': orange,
                  '': gray}
