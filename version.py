import git

GLOBAL_git_URL = "https://github.com/tamsen/FSAtoBoecheraMSI"
GLOBAL_public_version_info = [
    ["Differentiating between final calls and very confident calls, " +
     "so we only keep calls made with a very clear ladder. ", "Nov 14, 2022", "v1.0.0.0"],
    ["Adding background subtraction to peak finding for the ladder. " , "Nov 15, 2022", "v1.0.0.3"],
    ["Adding more smarts to ladder, so we have some expectation of where ladder should be, based on calibration data. ",
     "Nov 17, 2022", "v1.0.0.4"],
    ["Shifted A1 call by 0.5bp, made all the difference between paup and lxpxr call.", "Nov 22, 2022", "v1.0.0.5"],
    ["Outputs spreadsheet with species determinations", "Dec 01, 2022", "v1.0.0.6"],
    ["Extended to work with MW data, cleaned up input arguments", "Jan 10, 2023", "v1.0.1.0"],
    ["Extended to work with TM data", "Feb 15, 2023", "v1.0.2.0"],
    ["Tweeked ladder parameters to better deal with false peaks that came to close to real peaks on BF19",
     "Mar 20, 2023", "v1.0.3.0"],
    ["After that last change, more peaks would get called, causing more bugs, so I ended up undoing it",
     "Mar 24, 2023", "v1.0.3.1"],
    ["I noticed consistent cross contamination in PS1, between BF20 and A1 in TM's data," +
     " so I added a rule to remove it. Only affect's TM's data",
        "Feb 27, 2024", "v1.0.3.2"]
]



class version_info:

    def __init__(self):

        try:
            self.repo = git.Repo(search_parent_directories=True)
            self.sha = str(self.repo.head.commit)
            self.when = str(self.repo.head.commit.committed_datetime)

        except:
            self.repo = False
            self.sha = False
            self.when = False

        self.repo_url = GLOBAL_git_URL
        self.app_name = GLOBAL_git_URL.split("/")[-1]

        most_recent_update = GLOBAL_public_version_info[-1]
        self.version_num = most_recent_update[-1]
        self.public_comments = most_recent_update

    def most_recent_comment(self):
        return self.public_comments[-1][0]

    def to_string(self):

        data = ["\n\nGit repo:\t\t" + self.repo_url,
                "Version number:\t" + str(self.version_num),
                "Version notes:\t" + self.most_recent_comment()]

        if self.repo:
            data.append("Git changeset:\t" + self.sha + "\t" + self.when)

        return "\n".join(data)
