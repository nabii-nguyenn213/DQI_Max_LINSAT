from math import isqrt
from src.utils import read_csv
import numpy as np

def get_matrices(skillPath:str, slotPath:str, conflictPath:str, taskPath:str, quotaPath:str|None=None): 
    skills = read_csv(skillPath)[0]
    slots = read_csv(slotPath)[0]
    conflicts = read_csv(conflictPath)[0][:, 1:].astype(int)
    tasks = read_csv(taskPath)[0]
    subject_ids = read_csv(skillPath)[1][1:]
    slot_ids = read_csv(slotPath)[1][1:]
    if quotaPath is not None:
        quotas = read_csv(quotaPath)
        return skills, slots, conflicts, tasks, quotas, subject_ids, slot_ids
    return skills, slots, conflicts, tasks, subject_ids, slot_ids

def find_finite_field(lecturers:int # number of lecturers
                      ):
    q = max(2, lecturers)
    while True:
        if all(q%d!=0 for d in range(2, isqrt(q)+1)):
            return q 
        q+=1

def eligible_lecturers(skills:np.ndarray, slots:np.ndarray, tasks:np.ndarray, 
                       subject_ids, slot_ids, skillRequire=5, slotRequire=5):
    '''
    skills = [
            ['lecturer_1', 'skillsubject1_1', 'skillsubject2_1', ...], 
            ['lecturer_2', 'skillsubject1_2', 'skillsubject2_2', ...], 
            ...
            ['lecturer_L', 'skillsubject1_L', 'skillsubject2_L', ...]
            ] 

    slots = [
            ['lecturer_1', 'timeSlotID1_1', 'timeSlotID2_1', ...]
            ['lecturer_2', 'timeSlotID1_2', 'timeSlotID2_2', ...]
            ...
            ['lecturer_L', 'timeSlotID1_L', 'timeSlotID2_L', ...]
            ]

    tasks = [
            ['Class_1', 'Subject_1', 'timeSlotID_1', 'Room_1', 'Assigned Lecturer'], 
            ...
            ]
    '''
    skill_filter = np.asarray(skills[:, 1:], dtype=float) >= skillRequire
    slot_filter = np.asarray(slots[:, 1:], dtype=float) >= slotRequire

    subject_map = {name: i for i, name in enumerate(subject_ids)}
    slot_map = {name: i for i, name in enumerate(slot_ids)}

    lecturer_slot_map = {lecturer_id:i for i, lecturer_id in enumerate(slots[:, 0])}
    # print(lecturer_slot_map)

    eligible = dict()

    for task in tasks: 
        class_ = task[0]
        subject = task[1]
        time_slot = task[2]
        subject_idx = subject_map[subject]
        slot_idx = slot_map[time_slot]

        E = set()

        for i, lecturer_id in enumerate(skills[:, 0]):
            j = lecturer_slot_map[lecturer_id]
            qualified = skill_filter[i, subject_idx]
            available = slot_filter[j, slot_idx]
            if qualified and available: 
                E.add(i)

        eligible[class_] = E

    return eligible

def generate_conflicts(tasks: np.ndarray, slot_conflicts:np.ndarray, slot_ids): 
    '''
    slot_conflicts = [
                     [1, 0, 0, 1, ...], 
                     [0, 1, 0, 0, ...]
                     ... 
                     [0, 0, 0, 0, ...]
                     ]
        a square matrix nxn, n represent distinguished timetslotID, 
        where 0 if the slots is not conflict, 1 if the slots conflict
    '''
    n = len(tasks)
    k = len(slot_conflicts)
    slot_map = {str(slot_id): i for i, slot_id in enumerate(slot_ids)}
    # print("Slot map\n", slot_map)
    task_slots = [slot_map[str(task[2])] for task in tasks]
    # print("Task slot\n", task_slots) 
    conflicts = []
    for i in range(n): 
        for j in range(i+1, n): 
            a = task_slots[i]
            b = task_slots[j]
            # Same slot or overlapping slot
            if a == b or slot_conflicts[a, b]==1:
                conflicts.append((i, j))
    return conflicts

def formulate_max_linsat(q, eligible, conflicts): 
    n = len(eligible) # number of assigned lecturer 
    rows = []    # Constraint matrix
    allowed = [] # Acceptance sets
    weights = []
    labels  = []

    # Conflict constaints 
    for i, j in conflicts: 
        if not (0 <= i < n and 0 <= j < n and i!=j): 
            raise ValueError(f'Invalid conflict: {(i, j)}')
        row = np.zeros(n, dtype=int)
        row[i] = 1 
        row[j] = q - 1 # q-1 = -1 mod q
        rows.append(row)
        allowed.append([v!=0 for v in range(q)])
        weights.append(1.0)
        labels.append(f"conflict: {i}!={j}")

    # Eligibility constraints 
    for i, E in enumerate(eligible.values()): 
        E = set(E)
        if not E: 
            raise ValueError(f"Task {i} has no eligible lecturers")
        if len(E)==q: 
            continue 
        row = np.zeros(n, dtype=int)
        row[i]=1 
        rows.append(row)
        allowed.append([v in E for v in range(q)])
        weights.append(1.0)
        labels.append(f"eligible: {i}")
       
    # print("Rows =", rows)
    B = np.array(rows, dtype=int).reshape(-1, n)
    allowed = np.array(allowed, dtype=bool).reshape(-1, q)
    weights = np.array(weights, dtype=float)

    return B, allowed, weights, labels

def test_eligible(): 
    matrices = get_matrices(skillPath="datasets/tiny/InstructorSkill.csv", 
                    slotPath="datasets/tiny/InstructorSlot.csv", 
                    conflictPath="datasets/tiny/SlotConflict.csv", 
                    taskPath="datasets/tiny/Task.csv", 
                    quotaPath=None) 
    skills = matrices[0]
    slots = matrices[1]
    conflicts = matrices[2]
    tasks = matrices[3]
    subject_ids = matrices[4]
    slot_ids = matrices[5]
    print(eligible_lecturers(skills, slots, tasks, subject_ids, slot_ids))


def test_conflicts(): 
    matrices = get_matrices(skillPath="datasets/tiny/InstructorSkill.csv", 
                    slotPath="datasets/tiny/InstructorSlot.csv", 
                    conflictPath="datasets/tiny/SlotConflict.csv", 
                    taskPath="datasets/tiny/Task.csv", 
                    quotaPath=None) 
    skills = matrices[0]
    slots = matrices[1]
    conflicts = matrices[2]
    tasks = matrices[3]
    subject_ids = matrices[4]
    slot_ids = matrices[5]
    print("Conflicts:", generate_conflicts(tasks, conflicts, slot_ids))

if __name__ == "__main__": 
    # test_conflicts()
    matrices = get_matrices(skillPath="datasets/tiny/InstructorSkill.csv", 
                    slotPath="datasets/tiny/InstructorSlot.csv", 
                    conflictPath="datasets/tiny/SlotConflict.csv", 
                    taskPath="datasets/tiny/Task.csv", 
                    quotaPath=None) 
    skills = matrices[0]
    slots = matrices[1]
    conflicts = matrices[2]
    tasks = matrices[3]
    subject_ids = matrices[4]
    slot_ids = matrices[5]

    q = find_finite_field(skills.shape[0])
    eligible = eligible_lecturers(skills, slots, tasks, subject_ids, slot_ids)
    slot_conflicts = generate_conflicts(tasks, conflicts, slot_ids)

    B, allowed, weights, labels = formulate_max_linsat(q, eligible, slot_conflicts)

    print("Finite field:", q)
    print("Variables:", B.shape[1])
    print("Constraints:", B.shape[0])
    print("Constraint matrix\n", B)
    print("\nAcceptance Sets\n", allowed)
