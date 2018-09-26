        from matdb.utility import contract_absolute
        from matdb.io import save_dict_to_h5


#Temporary hack for line 361 in matdb/atoms.py
            if not path.isdir(data["folder"]):
                root = data["calc_contr_dir"]
                relpath = contract_absolute(data["folder"], root)
                parts = relpath.split('/')
                newname = parts[1] + '.' + parts[0].lower()
                parts[1] = newname
                newpath = '/'.join(parts)
                data["folder"] = path.join(root, newpath)
                print(relpath, data["folder"])
            #     with h5py.File(target,"w") as hf:
            #         save_dict_to_h5(hf,data,'/')
