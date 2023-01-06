if __name__=="__main__":
#     row = 3
#     column = 3
#     nums = row * column
#     num_list= []
#     for i in range(nums):
#         num_list.append(i)
#
#
#     class Solution(object):
#         def permute(self, nums):
#             result = []
#             self.permute_util(nums, 0, [], result)
#             return result
#
#         def permute_util(self, given_list, start, curr, result):
#             if start > len(given_list) - 1:
#                 # print(curr)
#                 result.append(curr)
#                 return
#             for i in range(start, len(given_list)):
#                 self.swap(given_list, start, start + (i - start))
#                 self.permute_util(given_list, start + 1, curr + [given_list[start]], result)
#                 # print(given_list)
#                 self.swap(given_list, start, start + (i - start))
#
#         def swap(self, nums, index1, index2):
#             temp = nums[index1]
#             nums[index1] = nums[index2]
#             nums[index2] = temp
#
#
#     ob1 = Solution()
#     counter = 0
#     permutations = (ob1.permute(num_list))
#     distinct = []
#     for i in permutations:
#         if distinct.count(i) == 0:
#             distinct.append(i)
#
#
#     for i in distinct:
#         if i[row*column-1] == 0:
#             counter += 1
#             print(i)
#     print(counter)
    nums = ["000000","110100", "011010", "001101" ]
    final =[]
    num = ''
    for i in nums:
        for j in nums:
            num = int(i)+int(j)
            num = str(num)
            new = ""
            if num not in final:
                final.append(num)
            print(num)
    print(final)