import git
from git import repo

GLOBAL_git_URL = "https://github.com/tamsen/FSAtoBoecheraMSI"
GLOBAL_version_num = "v1.0.0.0"
GLOBAL_public_comments = [
    ["Differentiating between final calls and very confident calls, " +
     "so we only keep calls made with a very clear ladder. ", "Nov 14, 2022"]
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
        self.version_num = GLOBAL_version_num
        self.public_comments = GLOBAL_public_comments

    def most_recent_comment(self):
        return self.public_comments[-1][0]

    def to_string(self):

        data = ["\n\nGit repo:\t\t" + self.repo_url,
                "Version number:\t" + str(self.version_num),
                "Version notes:\t" + self.most_recent_comment()]

        if (self.repo):
            data.append("Git changeset:\t" + self.sha + "\t" + self.when)

        return "\n".join(data)
