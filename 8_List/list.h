template <typename T, typename A = std::allocator<T>>
class List {
  private:
    struct BaseNode {
      BaseNode* next;
      BaseNode* prev;

      BaseNode() : next(this), prev(this) {}
    };

    struct Node: BaseNode {
      T value;
    };

    BaseNode fakeNode;
    size_t sz;

    using NodeAlloc = std::allocator_traits<A>::template rebind_alloc<Node>;
    using NodeTraits = std::allocator_traits<NodeAlloc>;
    [[no_unique_address]] NodeAlloc allocator;

    inline void push(BaseNode* node, const T& value) {
      Node* newNode = NodeTraits::allocate(allocator, 1);
      try {
        NodeTraits::construct(allocator, &newNode->value, value);
      } catch (...) {
        NodeTraits::deallocate(allocator, newNode, 1);
        throw;
      }
      ++sz;
      BaseNode* last = node->next;
      last->prev = node->next = newNode;
      newNode->prev = node;
      newNode->next = last;
    }

    inline void pop(BaseNode* bnode) {
      if (bnode == &fakeNode) return;
      Node* node = static_cast<Node*>(bnode);
      BaseNode* prv = node->prev;
      BaseNode* nxt = node->next;
      NodeTraits::destroy(allocator, &node->value);
      NodeTraits::deallocate(allocator, node, 1);
      --sz;
      prv->next = nxt;
      nxt->prev = prv;
    }

    inline void add(size_t n, const T& value) {
      for (size_t i = 0; i < n; ++i) {
        try {
          push_back(value);
        } catch (...) {
          while (!empty()) {
            pop_back();
          }
          throw;
        }
      }
    }

    void copy(const List& other) {
      int cnt = 0;
      for (Node* node = static_cast<Node*>(other.fakeNode.next); node != &other.fakeNode; node = static_cast<Node*>(node->next)) {
        try {
          push_back(node->value);
          ++cnt;
        } catch (...) {
          for (int i = 0; i < cnt; ++i) {
            pop_back();
          }
          throw;
        }
      }
    }

  public:
    using value_type = T;
    using allocator_type = A;

    void check() const {
      std::cout << "SIZE: " << sz << '\n';
      Node* node = static_cast<Node*>(fakeNode.next);
      for (size_t i = 0; i < sz; ++i, node = static_cast<Node*>(node->next)) {
        std::cout << "VALUE: " << node->value << " PREV: " << node->prev << " CUR: " << node << " NEXT: " << node->next << '\n';
      }
    }

    List() : fakeNode(), sz(0), allocator() {}

    List(size_t n) : List(n, T()) {}

    List(size_t n, const T& value) : List() { add(n, value); }

    List(const A& alctr) : fakeNode(), sz(0), allocator(alctr) {}

    List(size_t n, const A& alctr) : fakeNode(), sz(0), allocator(alctr) {
      for (size_t i = 0; i < n; ++i) {
        try {
          BaseNode* node = fakeNode.prev;
          Node* newNode = NodeTraits::allocate(allocator, 1);
          try {
            NodeTraits::construct(allocator, &newNode->value);
          } catch (...) {
            NodeTraits::deallocate(allocator, newNode, 1);
            throw;
          }
          ++sz;
          BaseNode* last = node->next;
          last->prev = node->next = newNode;
          newNode->prev = node;
          newNode->next = last;
        } catch (...) {
          while (!empty()) {
            pop_back();
          }
          throw;
        }
      }
    }

    List(size_t n, const T& value, const A& alctr) : fakeNode(), sz(0), allocator(alctr) { add(n, value); }

    List(const List& other) : fakeNode(), sz(0), allocator(NodeTraits::select_on_container_copy_construction(other.allocator)) {
      try {
        copy(other);
      } catch (...) {
        sz = 0;
      }
    }

    List& operator=(const List& other) {
      if (this == &other) return *this;
      size_t tmp = sz;
      try {
        copy(other);
        for (size_t i = 0; i < tmp; ++i) {
          pop_front();
        }
        if (NodeTraits::propagate_on_container_copy_assignment::value) {
          allocator = other.allocator;
        }
      } catch (...) {
        sz = tmp;
      }
      return *this;
    }

    ~List() noexcept {
      while (!empty()) {
        pop_back();
      }
    }

    void push_back(const T& value) { push(fakeNode.prev, value); }

    void pop_back() { pop(fakeNode.prev); }

    void push_front(const T& value) { push(&fakeNode, value); }

    void pop_front() { pop(fakeNode.next); }

    NodeAlloc get_allocator() const { return allocator; }

    size_t size() const { return sz; }

    bool empty() const { return sz == 0; }


    template<typename U, const int REV>
    class Iterator {
        friend class List;
      private:
        BaseNode* iter;

        explicit Iterator(BaseNode* iter) : iter(iter) {}

        explicit Iterator(const BaseNode* iter) : iter(const_cast<BaseNode*>(iter)) {}

      public:
        using value_type = U;
        using difference_type = int;
        using iterator_category = std::bidirectional_iterator_tag;
        using pointer = value_type*;
        using reference = value_type&;

        Iterator(const Iterator<U, REV>& other) : iter(other.iter) {}

        Iterator(const Iterator<U, -REV>& other) : iter(other.iter) {}

        operator Iterator<const U, REV>() {
          return Iterator<const U, REV>(iter);
        }

        Iterator& operator=(Iterator<U, REV> other) {
          iter = other.iter;
          return *this;
        }

        Iterator& operator++() {
          iter = (REV == 1) ? iter->next : iter->prev;
          return *this;
        }

        Iterator operator++(int) {
          Iterator temp(*this);
          ++*this;
          return temp;
        }

        Iterator& operator--() {
          iter = (REV == 1) ? iter->prev : iter->next;
          return *this;
        }

        Iterator operator--(int) {
          Iterator temp(*this);
          --*this;
          return temp;
        }

        Iterator<U, -REV> base() const {
          return Iterator<U, -REV>(iter->next);
        }

        bool operator==(const Iterator& other) const { return iter == other.iter; }

        bool operator!=(const Iterator& other) const { return !(*this == other); }

        U& operator*() { return static_cast<Node*>(iter)->value; }

        U* operator->() { return &(static_cast<Node*>(iter)->value); }

        const U& operator*() const { return static_cast<Node*>(iter)->value; }

        const U* operator->() const { return &(static_cast<Node*>(iter)->value); }
    };

    typedef Iterator<T, 1> iterator;
    typedef Iterator<const T, 1> const_iterator;
    typedef Iterator<T, -1> reverse_iterator;
    typedef Iterator<const T, -1> const_reverse_iterator;

    List<T, A>::iterator begin() { return iterator(fakeNode.next); }

    List<T, A>::iterator end() { return iterator(&fakeNode); }

    List<T, A>::const_iterator begin() const { return const_iterator(fakeNode.next); }

    List<T, A>::const_iterator end() const { return const_iterator(&fakeNode); }

    List<T, A>::reverse_iterator rbegin() { return reverse_iterator(fakeNode.prev); }

    List<T, A>::reverse_iterator rend() { return reverse_iterator(&fakeNode); }

    List<T, A>::const_reverse_iterator rbegin() const { return const_reverse_iterator(fakeNode.prev); }

    List<T, A>::const_reverse_iterator rend() const { return const_reverse_iterator(&fakeNode); }

    List<T, A>::const_iterator cbegin() const { return const_iterator(fakeNode.next); }

    List<T, A>::const_iterator cend() const { return const_iterator(&fakeNode); }

    List<T, A>::const_reverse_iterator crbegin() const { return const_reverse_iterator(fakeNode.prev); }

    List<T, A>::const_reverse_iterator crend() const { return const_reverse_iterator(&fakeNode); }

    template<typename U, const int REV>
    void insert(const Iterator<U, REV>& iter0, const T& value) { push(iter0.iter->prev, value); }

    template<typename U, const int REV>
    void erase(const Iterator<U, REV>& iter0) { pop(iter0.iter); }
};