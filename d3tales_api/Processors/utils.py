import numpy as np

MIN_SCAN_LEN = 5


# noinspection PyTypeChecker
def break_list(lst, return_brks=False):
    if len(lst) < 2:
        return lst
    direction = lst[1] - lst[0] / abs(lst[1] - lst[0])
    brks = []
    for i in range(1, len(lst)):
        previous_diff = lst[i] - lst[i - 1]
        if not previous_diff:
            continue
        i_direction = previous_diff / abs(previous_diff)
        if i_direction != direction:
            brks.append(i)
            direction = i_direction
    if return_brks:
        return brks
    return [lst[x:y] for x, y in zip([0] + brks, brks + [None])]


# noinspection PyTypeChecker
def break_scan_data(scan_data):
    scan_data = np.array(scan_data).tolist()
    volts = [v for v, i in scan_data]
    brks = break_list(volts, return_brks=True)
    return [scan_data[x:y] for x, y in zip([0] + brks, brks + [None])]


def mkdir_p(sftp, remote_directory):
    dir_path = str()
    for dir_folder in remote_directory.split("/"):
        if dir_folder == "":
            continue
        dir_path += r"/{0}".format(dir_folder)
        try:
            sftp.listdir(dir_path)
        except IOError:
            sftp.mkdir(dir_path)


def test_path(sftp, destination_path):
    folder_name = destination_path.split("/")[-1]
    folder_dir = "/".join(destination_path.split("/")[:-1])
    return True if folder_name in sftp.listdir(folder_dir) else False