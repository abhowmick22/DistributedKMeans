import difflib


def stringDistance(p1, p2):
    '''
    Takes two DNA strands and computes the distance between them.
    '''
    result = len(p1)         
    edits = 0
    difference = difflib.ndiff(p1, p2)
    for i,s in enumerate(difference):
        if s[0]=='-' or s[0]=='+':
          edits = edits + 1
    edits = edits/2
    return result - edits


print stringDistance("abcdefgh", "bbctefij")
