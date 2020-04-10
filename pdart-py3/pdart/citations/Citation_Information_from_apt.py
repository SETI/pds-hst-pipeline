import re
import xml.dom.minidom as md
from typing import List, Tuple, cast

# This is how we find a date embedded inside a comment
DATE_PATTERN1: re.Pattern = re.compile(
    r".*<!-.*Date.*[^0-9](19[7-9][0-9]|20[0-3][0-9])[^0-9].*"
)
DATE_PATTERN2: re.Pattern = re.compile(
    r"--.*Submission.*[^0-9](19[7-9][0-9]|20[0-3][0-9])[^0-9].*"
)

################################################################################


def Citation_Information_from_apt(
    filename: str,
) -> Tuple[int, str, int, List[str], str, str]:
    # Read file
    doc = md.parse(filename)

    # Get proposal number
    nodes: List[md.Node] = doc.getElementsByTagName("HSTProposal")
    propno = int(nodes[0].getAttribute("Phase2ID"))

    # Get category, cycle
    nodes = doc.getElementsByTagName("ProposalInformation")
    category: str = cast(str, nodes[0].getAttribute("Category"))
    cycle = int(nodes[0].getAttribute("Cycle"))

    # Get authors
    authors: List[str] = []
    for key in ("PrincipalInvestigator", "CoInvestigator"):
        nodes = doc.getElementsByTagName(key)
        for node in nodes:
            first = node.getAttribute("FirstName")
            last = node.getAttribute("LastName")
            middle = node.getAttribute("MiddleInitial")

            if middle:
                author = first + " " + middle + " " + last
            else:
                author = first + " " + last

            authors.append(author)

    # Get title
    nodes = doc.getElementsByTagName("Title")
    title: str = cast(md.Text, nodes[0].childNodes[0]).data

    words = title.split()  # clean up whitespace
    title = " ".join(words)

    # Try to get the year from a "Date" comment or Submission line
    with open(filename) as f:
        recs = f.readlines()

    year = 0
    for rec in recs:
        if year != 0:
            break

        for pattern in (DATE_PATTERN1, DATE_PATTERN2):
            match = re.match(pattern, rec)
            if match:
                year = int(match.group(1))
                break

    # See if there is a later year in a program constraint
    for key in ("Start", "End"):
        nodes = doc.getElementsByTagName(key)
        for node in nodes:
            year_str = node.getAttribute("Year")
            if year_str:
                year = max(year, int(year_str))

    if year == 0:
        raise ValueError("missing year in " + filename)

    return (propno, category, cycle, authors, title, str(year))
