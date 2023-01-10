
def show_example_usage():
    print("o=/home/tamsen/Data/Eton i=./data/FSAlist.txt panel=./data/Panel_mw.xml truth=./data/TruthData.xml")

def do_parsing(arg_list):


    truth_file = False
    upgraded_args= {}
    rules="Tamsen"

    for arg in arg_list:
        if "=" in arg:
            splat=arg.split("=")
            key=splat[0]
            value=splat[1]
            upgraded_args[key]=value

    if "truth" in upgraded_args:
        truth_file = upgraded_args["truth"]

    if "rules" in upgraded_args:
        rules= upgraded_args["rules"]

    if rules != "Tamsen":
        if rules != "MW":
            print("Sorry, you need to use MW or Tamsen's rules. Try again.")
            raise NotImplementedError()

    return [upgraded_args["o"],upgraded_args["i"],upgraded_args["panel"],truth_file, rules]